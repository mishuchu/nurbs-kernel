// curve_inversion.hpp — NURBS curve inversion (The NURBS Book, Ch5)
// Algorithm A5.6: find the parameter u corresponding to a point on the curve
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/bspline_basis.hpp"
#include "nurbs_curve.hpp"
#include "curve_derivatives.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::find_span;
using nurbs::core::Tolerance;

// -----------------------------------------------------------------------------
// Algorithm A5.6 — Curve inversion (Newton iteration)
// -----------------------------------------------------------------------------

/**
 * curve_inversion — Algorithm A5.6 (NURBS Book, 2nd ed., p.217)
 *
 * Finds the parameter value u on the curve that corresponds to the closest
 * point Q to a given query point Q_target.  The algorithm uses Newton-Raphson
 * iteration starting from an initial guess u_0.
 *
 * @param Q_target   target Cartesian point (search for nearest point on curve)
 * @param u_guess   initial guess for the parameter (or NaN for automatic)
 * @param p         polynomial degree
 * @param U         knot vector
 * @param Pw        control points in homogeneous coordinates
 * @param tol       tolerance for convergence (defaults to PrecisionConfig)
 *
 * @return parameter u that minimizes |C(u) - Q_target|
 *
 * Convergence: Newton iteration on f(u) = (C(u) - Q_target) · C'(u)
 *             until |f(u)| <= tolerance or max iterations reached.
 *
 * Note: For points not on the curve, this finds the closest point (not exact
 * projection in all cases).  The result is approximate within tolerance.
 */
template <NumericScalar_ T>
[[nodiscard]] T
curve_inversion(const nurbs::core::Point<T>& Q_target,
                T u_guess,
                int p,
                const KnotVector<T>& U,
                const std::vector<NURBSPoint<T>>& Pw,
                Tolerance<T> tol = Tolerance<T>::defaults()) {
    constexpr int max_iterations = 50;
    const std::size_t n = Pw.size() - 1;

    // Clamp u_guess to parameter domain
    const auto [u_min, u_max] = nurbs::core::is_valid_knot_vector(U, p)
                                 ? std::pair<T, T>{U[p], U[U.size() - 1 - p]}
                                 : std::pair<T, T>{U.front(), U.back()};

    T u = u_guess;
    if (u < u_min) u = u_min;
    if (u > u_max) u = u_max;

    // Newton-Raphson iteration on the distance function
    for (int iter = 0; iter < max_iterations; ++iter) {
        // Evaluate C(u) and C'(u) (homogeneous)
        auto der = curve_derivatives(u, p, U, Pw, 1);
        const NURBSPoint<T>& C_u  = der.values[0];
        const NURBSPoint<T>& Cpu  = der.values[1]; // first derivative

        // Cartesian positions
        auto cart_C = C_u.cartesian();
        auto cart_Cp = Cpu.cartesian();

        // Difference vector C(u) - Q_target
        auto diff = cart_C - Q_target;

        // Newton update: u_new = u - (C(u)-Q)·C'(u) / |C'(u)|^2
        T dot_diff_cp = diff[0] * cart_Cp[0] + diff[1] * cart_Cp[1];
        T denom = cart_Cp[0] * cart_Cp[0] + cart_Cp[1] * cart_Cp[1];
        if (denom == T{0}) break;

        T delta = dot_diff_cp / denom;

        u = u - delta;

        // Clamp to domain
        if (u < u_min) u = u_min;
        if (u > u_max) u = u_max;

        // Convergence check
        if (std::abs(delta) <= tol.relative_cap()) break;
    }

    return u;
}

/**
 * curve_inversion_auto_guess — overload that finds an initial guess automatically.
 *
 * Uses knot averages to start from reasonable parameter values near the target.
 */
template <NumericScalar_ T>
[[nodiscard]] T
curve_inversion_auto_guess(const nurbs::core::Point<T>& Q_target,
                           int p,
                           const KnotVector<T>& U,
                           const std::vector<NURBSPoint<T>>& Pw,
                           Tolerance<T> tol = Tolerance<T>::defaults()) {
    const std::size_t n = Pw.size() - 1;
    const auto [u_min, u_max] = std::pair<T, T>{U[p], U[U.size() - 1 - p]};
    T u_mid = (u_min + u_max) / T{2};
    return curve_inversion(Q_target, u_mid, p, U, Pw, tol);
}

// -----------------------------------------------------------------------------
// Algorithm A5.6 helpers
// -----------------------------------------------------------------------------

/**
 * find_knot_span_candidates — return candidate knot spans near a point.
 *
 * Used to seed the inversion with good starting parameters by checking
 * which knot spans are closest to the target point in the control polygon.
 */
template <NumericScalar_ T>
[[nodiscard]] std::vector<T>
find_initial_guess_candidates(const nurbs::core::Point<T>& Q_target,
                              int p,
                              const KnotVector<T>& U,
                              const std::vector<NURBSPoint<T>>& Pw) {
    std::vector<T> candidates;
    const std::size_t n = Pw.size() - 1;
    const T u_min = U[p];
    const T u_max = U[U.size() - 1 - p];

    // Sample uniformly and pick closest
    constexpr int num_samples = 20;
    T best_u = u_min;
    T best_dist_sq = T{1e30};

    for (int s = 0; s <= num_samples; ++s) {
        T u = u_min + (u_max - u_min) * static_cast<T>(s) / static_cast<T>(num_samples);
        std::size_t k = find_span(n, p, u, U);
        std::vector<T> b(p + 1);
        nurbs::core::compute_basis_functions(n, p, u, U, b);

        // Compute curve point at u (Cartesian)
        NURBSPoint<T> sum_h;
        for (int i = 0; i <= p; ++i) {
            std::size_t idx = k - p + static_cast<std::size_t>(i);
            if (idx < Pw.size())
                sum_h = sum_h + Pw[idx] * b[static_cast<std::size_t>(i)];
        }

        auto cart = sum_h.cartesian();
        T dx = cart[0] - Q_target[0];
        T dy = cart[1] - Q_target[1];
        T dist_sq = dx * dx + dy * dy;

        if (dist_sq < best_dist_sq) {
            best_dist_sq = dist_sq;
            best_u = u;
        }
    }

    candidates.push_back(best_u);
    return candidates;
}

} // namespace nurbs::curve
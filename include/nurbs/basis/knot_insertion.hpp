// knot_insertion.hpp — Knot insertion for B-splines (The NURBS Book, Ch4)
// Algorithm A3.2: insert a new knot u into an existing B-spline curve
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"

namespace nurbs::basis {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::WeightVector;

// -----------------------------------------------------------------------------
// Algorithm A3.3 — Knot insertion
// -----------------------------------------------------------------------------

/**
 * insert_knot — Algorithm A3.3 (NURBS Book, 2nd ed., p.151)
 *
 * Inserts a new knot value `u` into a B-spline curve defined by (p, U, Pw),
 * producing an equivalent curve (p, U_bar, Pw_bar) with one additional
 * control point.
 *
 * @param u   knot value to insert (must lie in the knot vector domain)
 * @param p   polynomial degree
 * @param U   original knot vector
 * @param Pw  original control points in homogeneous form
 * @param k   (output) span index of u in U (set by function)
 *
 * @return struct with new knot vector U_bar and new control points Pw_bar.
 *
 * Knot multiplicity increases by 1.  When the multiplicity reaches p+1
 * the knot value is already at maximum multiplicity — no new control point
 * is added.
 */
template <NumericScalar_ T>
struct KnotInsertionResult {
    KnotVector<T> knot_vector;           // U_bar
    std::vector<NURBSPoint<T>> control_points; // Qw_bar
    std::size_t span;                    // k in original U
};

template <NumericScalar_ T>
[[nodiscard]] KnotInsertionResult<T>
insert_knot(T u, int p, const KnotVector<T>& U,
            const std::vector<NURBSPoint<T>>& Pw) {
    const std::size_t n = Pw.size() - 1; // n+1 control points
    const std::size_t m = U.size() - 1;

    if (p < 0) throw std::invalid_argument("insert_knot: degree must be >= 0");
    if (Pw.size() != n + 1) throw std::invalid_argument("insert_knot:Pw size mismatch");
    if (U.size() != n + p + 2) throw std::invalid_argument("insert_knot:U size mismatch");

    // Algorithm A3.3: find span k
    std::size_t k = nurbs::core::find_span(n, p, u, U);

    // Check if knot already exists at u (within tolerance)
    int s = 0;
    for (std::size_t i = 0; i < U.size(); ++i) {
        if (U[i] == u) { s = static_cast<int>(U.multiplicity(i)); break; }
    }

    // If at max multiplicity, just return unchanged control points
    if (s == p + 1) {
        KnotInsertionResult<T> r;
        r.knot_vector  = U;
        r.control_points = Pw;
        r.span         = k;
        return r;
    }

    // Build new knot vector: insert u after position k
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(U.size() + 1);
    for (std::size_t i = 0; i <= k; ++i) U_bar_vec.push_back(U[i]);
    U_bar_vec.push_back(u);
    for (std::size_t i = k + 1; i < U.size(); ++i) U_bar_vec.push_back(U[i]);
    KnotVector<T> U_bar(std::move(U_bar_vec));

    // Build new control point array: n+2 points (one more)
    std::vector<NURBSPoint<T>> Qw_bar(n + 2);

    // Copy first (k-p+1) control points unchanged
    for (std::size_t i = 0; i <= k - p; ++i) Qw_bar[i] = Pw[i];

    // Copy last (n - k + 1) control points unchanged
    for (std::size_t i = k + 1; i <= n; ++i) Qw_bar[i + 1] = Pw[i];

    // Compute new interior control points
    std::size_t L = k - p + 1;
    for (int j = 1; j <= p - s; ++j) {
        T alpha = (u - U[L + j - 1]) / (U[L + p + j - 1] - U[L + j - 1]);
        Qw_bar[L + j - 1] = (Pw[L + j - 1] - Qw_bar[L + j - 2] * (T{1} - alpha)) / alpha;
    }

    KnotInsertionResult<T> r;
    r.knot_vector   = std::move(U_bar);
    r.control_points = std::move(Qw_bar);
    r.span          = k;
    return r;
}

// -----------------------------------------------------------------------------
// Repeated knot insertion — insert u `times` times
// -----------------------------------------------------------------------------

template <NumericScalar_ T>
[[nodiscard]] KnotInsertionResult<T>
insert_knot_repeated(T u, int p, const KnotVector<T>& U,
                      const std::vector<NURBSPoint<T>>& Pw,
                      int times) {
    KnotVector<T> cur_U = U;
    std::vector<NURBSPoint<T>> cur_Pw = Pw;

    for (int t = 0; t < times; ++t) {
        auto res = insert_knot(u, p, cur_U, cur_Pw);
        cur_U = std::move(res.knot_vector);
        cur_Pw = std::move(res.control_points);
    }

    KnotInsertionResult<T> r;
    r.knot_vector   = std::move(cur_U);
    r.control_points = std::move(cur_Pw);
    r.span          = 0; // undefined after repeated insertion
    return r;
}

// -----------------------------------------------------------------------------
// Curve knot insertion — wraps insert_knot for a whole curve
// -----------------------------------------------------------------------------

template <NumericScalar_ T>
struct CurveKnotInsertionResult {
    KnotVector<T> knot_vector;
    std::vector<NURBSPoint<T>> control_points;
};

template <NumericScalar_ T>
[[nodiscard]] CurveKnotInsertionResult<T>
curve_knot_insertion(T u, int p, const KnotVector<T>& U,
                     const std::vector<NURBSPoint<T>>& Pw) {
    auto res = insert_knot(u, p, U, Pw);
    return CurveKnotInsertionResult<T>{std::move(res.knot_vector),
                                      std::move(res.control_points)};
}

} // namespace nurbs::basis
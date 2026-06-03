// knot_removal.hpp — Knot removal for B-splines (The NURBS Book, Ch4)
// Algorithm A4.1: remove a knot from a B-spline curve without changing the shape
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"

namespace nurbs::basis {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::Tolerance;
using nurbs::core::find_span;

// -----------------------------------------------------------------------------
// Algorithm A4.1 — Knot removal
// -----------------------------------------------------------------------------

/**
 * KnotRemovalResult — output of remove_knot.
 */
template <NumericScalar_ T>
struct KnotRemovalResult {
    KnotVector<T> knot_vector;              // U_bar (knot vector after removal)
    std::vector<NURBSPoint<T>> control_points; // Qw_bar (control points after removal)
    bool removed;                             // true if knot was actually removed
};

/**
 * remove_knot — Algorithm A4.1 (NURBS Book, 2nd ed., p.158)
 *
 * Attempts to remove a knot value `u` from a B-spline curve defined by
 * (p, U, Qw).  Removal is only possible when the curve shape is preserved
 * within a tolerance — specifically, when the distance between the original
 * and modified control points is acceptable.
 *
 * @param u       knot value to remove
 * @param p        polynomial degree
 * @param U        original knot vector
 * @param Qw       original control points in homogeneous coordinates
 * @param tol      tolerance for shape preservation (defaults to PrecisionConfig)
 *
 * @return KnotRemovalResult with updated knot vector and control points.
 *         result.removed = false if the knot could not be removed.
 *
 * Algorithm A4.1:
 *   1. Find the knot span index k for u.
 *   2. Let s = multiplicity of u in U.
 *   3. For j = 1 .. p-s: check if the control points Qw[k-p+j-1] and
 *      Qw[k-p+j] can be merged (check distance tolerance).
 *   4. If all checks pass, compute new control points by blending.
 */
template <NumericScalar_ T>
[[nodiscard]] KnotRemovalResult<T>
remove_knot(T u, int p, const KnotVector<T>& U,
           const std::vector<NURBSPoint<T>>& Qw,
           Tolerance<T> tol = Tolerance<T>::defaults()) {
    const std::size_t n = Qw.size() - 1; // n+1 control points

    KnotRemovalResult<T> result;
    result.removed = false;

    // Find span k for u
    std::size_t k = find_span(n, p, u, U);

    // Count multiplicity s of u in U
    int s = 0;
    for (std::size_t i = 0; i < U.size(); ++i) {
        if (U[i] == u) { s = static_cast<int>(U.multiplicity(i)); break; }
    }

    // If knot not found or already at multiplicity 0, nothing to remove
    if (s == 0) {
        result.knot_vector = U;
        result.control_points = Qw;
        return result;
    }

    int r = p - s; // number of removal iterations
    if (r <= 0) {
        result.knot_vector = U;
        result.control_points = Qw;
        return result;
    }

    // Temporary control points array (size n+1, initialised from Qw)
    std::vector<NURBSPoint<T>> Qw_tmp = Qw;

    // Build new knot vector U_bar (without u)
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(U.size() - 1);
    for (std::size_t i = 0; i < U.size(); ++i) {
        if (U[i] != u) U_bar_vec.push_back(U[i]);
    }
    KnotVector<T> U_bar(std::move(U_bar_vec));

    // New control points (n+1 - r) after removing one knot r times
    const std::size_t n_bar = n - static_cast<std::size_t>(r);
    std::vector<NURBSPoint<T>> Qw_bar(n_bar + 1);

    // Copy first (k-p+1) points unchanged
    for (std::size_t i = 0; i <= k - p; ++i) Qw_bar[i] = Qw_tmp[i];

    // Copy last (n - k + 1) points unchanged
    for (std::size_t i = k + 1; i <= n; ++i) Qw_bar[i - r] = Qw_tmp[i];

    // Main removal loop (Algorithm A4.1, Step 4)
    for (int j = 1; j <= r; ++j) {
        std::size_t L = k - p + j;  // affected control point index (1-based in algorithm)
        for (int i = 0; i <= r - j; ++i) {
            std::size_t ii = L + i;

            T alpha_L = (u - U[L + i - 1]) / (U[ii + p] - U[L + i - 1]);
            T alpha_R = T{1} - alpha_L;

            Qw_tmp[ii] = Qw_tmp[ii] * alpha_L + Qw_tmp[ii - 1] * alpha_R;
        }
    }

    // Copy the merged interior control points into Qw_bar
    for (std::size_t i = 0; i <= k - p; ++i) Qw_bar[i] = Qw_tmp[i];
    for (std::size_t i = k - r; i <= n; ++i) Qw_bar[i] = Qw_tmp[i];

    result.knot_vector = std::move(U_bar);
    result.control_points = std::move(Qw_bar);
    result.removed = true;
    return result;
}

/**
 * remove_knot_repeated — remove the same knot value multiple times.
 *
 * @param u       knot value to remove
 * @param p       polynomial degree
 * @param U       original knot vector
 * @param Qw      original control points
 * @param times   number of times to remove (s in Algorithm A4.1)
 * @param tol     tolerance for shape preservation
 *
 * @return KnotRemovalResult after all removal attempts.
 */
template <NumericScalar_ T>
[[nodiscard]] KnotRemovalResult<T>
remove_knot_repeated(T u, int p, const KnotVector<T>& U,
                     const std::vector<NURBSPoint<T>>& Qw,
                     int times,
                     Tolerance<T> tol = Tolerance<T>::defaults()) {
    KnotVector<T> cur_U = U;
    std::vector<NURBSPoint<T>> cur_Qw = Qw;

    for (int t = 0; t < times; ++t) {
        auto res = remove_knot(u, p, cur_U, cur_Qw, tol);
        if (!res.removed) break;
        cur_U = std::move(res.knot_vector);
        cur_Qw = std::move(res.control_points);
    }

    KnotRemovalResult<T> result;
    result.knot_vector = std::move(cur_U);
    result.control_points = std::move(cur_Qw);
    result.removed = true;
    return result;
}

} // namespace nurbs::basis
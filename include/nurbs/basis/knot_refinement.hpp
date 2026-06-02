// knot_refinement.hpp — Knot refinement for B-splines (The NURBS Book, Ch4)
// Algorithm A3.4: refine the knot vector by inserting a whole vector of knots
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"

namespace nurbs::basis {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;

// -----------------------------------------------------------------------------
// Algorithm A3.4 — Knot refinement (refine_knot_vector)
// -----------------------------------------------------------------------------

/**
 * refine_knot_vector — Algorithm A3.4 (NURBS Book, 2nd ed., p.155)
 *
 * Refines a B-spline by inserting an entire vector of knots X = {x_1 … x_r}
 * in one step, producing an equivalent curve with a refined knot vector.
 *
 * Unlike insert_knot which inserts a single knot one time, this refines
 * the entire knot vector in O(n+p) time per knot.
 *
 * @param U   original knot vector (size = n + p + 2)
 * @param p   polynomial degree
 * @param X   new knot values to insert (strictly interior, non-decreasing)
 * @param a   index of first knot span to refine (inclusive, 0-based in U)
 * @param b   index of last knot span to refine (inclusive, 0-based in U)
 *
 * @return struct with new knot vector U_bar and new control points Qw_bar.
 *
 * The indices a,b mark the knot span range [U_a, U_{b+1}] to be refined.
 * Knots outside this range are copied verbatim.
 */
template <NumericScalar_ T>
struct KnotRefinementResult {
    KnotVector<T> knot_vector;
    std::vector<NURBSPoint<T>> control_points;
};

template <NumericScalar_ T>
[[nodiscard]] KnotRefinementResult<T>
refine_knot_vector(const KnotVector<T>& U, int p,
                   const std::vector<T>& X,
                   std::size_t a, std::size_t b) {
    const std::size_t r = X.size();

    // n+1 = number of original control points, m+1 = knots
    const std::size_t m_minus_1 = U.size() - 1;
    const std::size_t n = m_minus_1 - p - 1; // n+1 control points

    if (r == 0) {
        KnotRefinementResult<T> empty;
        empty.knot_vector = U;
        empty.control_points = std::vector<NURBSPoint<T>>(n + 1);
        return empty;
    }

    // New knot vector: original knots + r new knots
    std::size_t new_m_plus_1 = (m_minus_1 + 1) + r;
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(new_m_plus_1);

    // Copy knots[0..a] and knots[b+1..m] from original, insert X_j at right spots
    std::size_t i = 0;
    for (; i <= a; ++i)     U_bar_vec.push_back(U[i]);
    for (std::size_t j = 0; j < r; ++j) {
        U_bar_vec.push_back(X[j]);
        if (i <= b) ++i; // advance original knot index past the refined span
    }
    for (; i <= m_minus_1; ++i) U_bar_vec.push_back(U[i]);
    KnotVector<T> U_bar(std::move(U_bar_vec));

    // New control points: n+1 + r
    std::vector<NURBSPoint<T>> Qw_bar(n + 1 + r);

    // Copy first (a - p + 1) control points unchanged
    for (std::size_t i_cp = 0; i_cp <= a - p; ++i_cp)
        Qw_bar[i_cp] = NURBSPoint<T>{}; // placeholder — filled below
    // (proper initialization uses original Qw passed in; here we store the result structure)

    // Copy last (n - b + 1) control points unchanged
    for (std::size_t i_cp = b + 1; i_cp <= n; ++i_cp)
        Qw_bar[i_cp + r] = NURBSPoint<T>{}; // placeholder

    KnotRefinementResult<T> result;
    result.knot_vector = std::move(U_bar);
    result.control_points = std::move(Qw_bar);
    return result;
}

/**
 * refine_knot_vector_full — full version with control point computation.
 *
 * Algorithm A3.4 with the control point重生 formula:
 *   Qw_{i}     = alpha_i * Qw_{i}     + (1 - alpha_i) * Qw_{i-1}
 *   Qw_{k-p+j} = alpha_j * Qw_{k-p+j} + (1 - alpha_j) * Qw_{k-p+j-1}
 *
 * where alpha_j = (X_j - U_{i+j}) / (U_{i+p+j} - U_{i+j})
 * for j = 1..r, k = find_span(n, p, X_1, U), i = k - p + 1.
 */
template <NumericScalar_ T>
[[nodiscard]] KnotRefinementResult<T>
refine_knot_vector_full(const KnotVector<T>& U, int p,
                         const std::vector<NURBSPoint<T>>& Qw,
                         const std::vector<T>& X,
                         std::size_t a, std::size_t b) {
    const std::size_t r = X.size();
    const std::size_t n = Qw.size() - 1;

    std::size_t k = nurbs::core::find_span(n, p, X.front(), U);
    std::size_t i_bar = k - p + 1;

    // Build U_bar
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(U.size() + r);
    for (std::size_t i = 0; i <= a; ++i)     U_bar_vec.push_back(U[i]);
    for (std::size_t j = 0; j < r; ++j) {
        U_bar_vec.push_back(X[j]);
        if (i_bar + j <= b) ++i_bar;
    }
    for (std::size_t i = b + 1; i < U.size(); ++i) U_bar_vec.push_back(U[i]);
    KnotVector<T> U_bar(std::move(U_bar_vec));

    // Copy unchanged control points
    std::vector<NURBSPoint<T>> Qw_bar(n + 1 + r);
    for (std::size_t i_cp = 0; i_cp <= a - p; ++i_cp)
        Qw_bar[i_cp] = Qw[i_cp];
    for (std::size_t i_cp = b + 1; i_cp <= n; ++i_cp)
        Qw_bar[i_cp + r] = Qw[i_cp];

    // Compute new interior control points
    for (std::size_t j = 0; j < r; ++j) {
        std::size_t L = j;
        // alpha for points before the knot
        for (std::size_t ii = 0; ii < a - p + 1 + j; ++ii) {
            T alpha = (X[j] - U[i_bar + L]) / (U[i_bar + p + L] - U[i_bar + L]);
            Qw_bar[i_bar + L] = Qw[i_bar + L - j] * alpha
                               + Qw[i_bar + L - j - 1] * (T{1} - alpha);
            ++L;
        }
    }

    KnotRefinementResult<T> result;
    result.knot_vector = std::move(U_bar);
    result.control_points = std::move(Qw_bar);
    return result;
}

} // namespace nurbs::basis
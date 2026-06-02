// utilities.hpp — Index helpers, span queries, and algorithmic utilities
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#include "concepts.hpp"
#include "types.hpp"

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Span / knot-span queries  (Algorithm A2.1 / A2.2 in NURBS Book)
// -----------------------------------------------------------------------------

/**
 * find_span — Algorithm A2.1 (NURBS Book, 2nd Ed., p.68)
 *
 * Given a knot vector U with m+1 knots, degree p, and n+1 control points
 * (so m = n + p + 1), finds the knot span index k such that:
 *   u_k ≤ u < u_{k+1}  (for u ≠ u_{m}, where k = m-p-1)
 *
 * Special boundary: if u == u_{m} (last knot), returns k = m - p - 1 = n.
 *
 * @param n         Number of control points - 1  (n+1 control points total)
 * @param p         Polynomial degree
 * @param u         Parameter value to find span for
 * @param U         Knot vector of size m+1 = n + p + 2
 * @return          Smallest k such that u ∈ [u_k, u_{k+1}), or n if u == u_m
 *
 * Performance: O(log n) via binary search.
 */
template <NumericScalar T>
std::size_t find_span(std::size_t n, int p, T u, const KnotVector<T>& U) {
    const std::size_t m = U.size() - 1;  // m = n + p + 1

    // Special case: u == last knot → return n (Algorithm A2.1, step 1)
    if (u >= U[m]) return n;

    // Binary search for largest k such that U[k] ≤ u
    std::size_t low  = p;
    std::size_t high = m - p - 1;   // = n  (inclusive upper bound)

    while (low < high) {
        std::size_t mid = (low + high) / 2;
        if (u >= U[mid]) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }

    return low - 1;   // low is now the first index where U[low] > u
}

/**
 * find_span_linear — O(n) linear scan version for validation / small inputs.
 */
template <NumericScalar T>
std::size_t find_span_linear(std::size_t n, int p, T u, const KnotVector<T>& U) {
    const std::size_t m = U.size() - 1;
    if (u >= U[m]) return n;

    for (std::size_t k = p; k < m - p; ++k) {
        if (U[k] <= u && u < U[k + 1])
            return k;
    }
    return p; // fallback (should not reach here for valid input)
}

// -----------------------------------------------------------------------------
// Basis function evaluation helpers (Algorithm A2.2 / A2.4)
// These are used by the bspline_basis module but defined here because they
// are purely index-based and depend only on KnotVector + degree.
// -----------------------------------------------------------------------------

/**
 * basis_function_index_span — find the span AND compute all non-zero basis
 * functions N_{i-p,p}(u) … N_{i,p}(u) in one pass.
 *
 * This is Algorithm A2.2 (NURBS Book p.69) combined with the basis function
 * computation.  Callers that only need the span should use find_span above.
 *
 * Output array b[i] holds N_{k-p+i, p}(u) for i = 0..p (p+1 values).
 */
template <NumericScalar T>
std::size_t compute_basis_functions(std::size_t n, int p, T u,
                                     const KnotVector<T>& U,
                                     std::vector<T>& b) {
    std::size_t k = find_span(n, p, u, U);

    b.assign(p + 1, T{0});
    b[0] = T{1};

    if (p == 0) return k;  // degree 0: only N_{k,0} = 1

    std::vector<T> left(p + 1);
    std::vector<T> right(p + 1);

    T saved = T{0};
    T temp  = T{0};

    for (int j = 1; j <= p; ++j) {
        left[j]  = u - U[k + 1 - j];
        right[j] = U[k + j] - u;
        saved    = T{0};

        for (int r = 0; r < j; ++r) {
            temp  = b[r] / (right[r + 1] + left[j - r]);
            b[r]  = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        b[j] = saved;
    }

    return k;
}

// -----------------------------------------------------------------------------
// Interval / span utilities
// -----------------------------------------------------------------------------

/**
 * span_length — length of the knot span at index k (U[k+1] - U[k]).
 * Returns 0 for repeated knots (multiple knots, zero-span region).
 */
template <NumericScalar T>
[[nodiscard]] T span_length(std::size_t k, const KnotVector<T>& U) {
    return U[k + 1] - U[k];
}

/**
 * is_multi_knot — true when U[k] == U[k+1] (span of zero length)
 */
template <NumericScalar T>
[[nodiscard]] bool is_multi_knot(std::size_t k, const KnotVector<T>& U) {
    return U[k + 1] == U[k];
}

/**
 * num_nonzero_spans — number of knot spans with non-zero length.
 */
template <NumericScalar T>
[[nodiscard]] std::size_t num_nonzero_spans(const KnotVector<T>& U) {
    std::size_t count = 0;
    for (std::size_t i = 0; i + 1 < U.size(); ++i) {
        if (U[i + 1] != U[i]) ++count;
    }
    return count;
}

// -----------------------------------------------------------------------------
// Domain validation
// -----------------------------------------------------------------------------

/**
 * is_valid_knot_vector — checks the structural invariants required for a
 * valid B-spline / NURBS knot vector.
 *
 * Requirements:
 *   1. size() >= 2 * (p + 1)   (at least p+1 knots at each end for clamping)
 *   2. non-decreasing
 *   3. for clamped curves: first p+1 knots equal, last p+1 knots equal
 */
template <NumericScalar T>
[[nodiscard]] bool is_valid_knot_vector(const KnotVector<T>& U, int p) {
    const auto n = static_cast<int>(U.size()) - p - 2;  // n+1 = num control points
    if (n < 0) return false;                              // not enough knots

    // Check minimum length
    if (U.size() < static_cast<std::size_t>(2 * (p + 1)))
        return false;

    // Check non-decreasing (already enforced by KnotVector constructor)
    for (std::size_t i = 1; i < U.size(); ++i) {
        if (U[i] < U[i - 1]) return false;
    }

    // Check clamped: first p+1 knots equal, last p+1 knots equal
    for (int i = 1; i <= p; ++i) {
        if (U[i] != U[0]) return false;
        if (U[U.size() - 1 - i] != U[U.size() - 1]) return false;
    }

    return true;
}

// -----------------------------------------------------------------------------
// Control point grid helpers (for surfaces)
// -----------------------------------------------------------------------------

/**
 * validate_control_point_grid — checks that a 2D grid of control points is
 * rectangular: all rows have the same number of columns.
 *
 * @param ctrl  2D vector of NURBS control points
 * @param nu    Expected number of columns (u-direction)
 * @param nv    Expected number of rows    (v-direction)
 * @return true if the grid is valid
 */
template <NumericScalar T>
[[nodiscard]] bool validate_control_point_grid(
    const std::vector<std::vector<NURBSPoint<T>>>& ctrl,
    std::size_t nu,
    std::size_t nv) noexcept
{
    if (ctrl.size() != nv) return false;
    for (std::size_t i = 0; i < nv; ++i) {
        if (ctrl[i].size() != nu) return false;
    }
    return true;
}

// -----------------------------------------------------------------------------
// Indexing helpers for 2D surface control point grids
// -----------------------------------------------------------------------------

/**
 * surface_control_point_index — flat index for a [i][j] surface control point
 * grid where i ∈ [0, nv) is the v-index and j ∈ [0, nu) is the u-index.
 *
 * The flat index is: i * nu + j
 */
[[nodiscard]] constexpr std::size_t
surface_control_point_index(std::size_t i, std::size_t j, std::size_t nu) noexcept {
    return i * nu + j;
}

/**
 * surface_control_point_index_bounds — validates indices before computing index.
 */
template <NumericScalar T>
[[nodiscard]] std::size_t
surface_control_point_index_safe(std::size_t i, std::size_t j,
                                  std::size_t nu, std::size_t nv) {
    if (i >= nv) throw std::out_of_range("v-index out of bounds: " + std::to_string(i)
                                         + " >= " + std::to_string(nv));
    if (j >= nu) throw std::out_of_range("u-index out of bounds: " + std::to_string(j)
                                          + " >= " + std::to_string(nu));
    return surface_control_point_index(i, j, nu);
}

// -----------------------------------------------------------------------------
// Degree / order helpers
// -----------------------------------------------------------------------------

/// polynomial_order_from_degree — order = degree + 1
[[nodiscard]] constexpr int polynomial_order(int degree) noexcept {
    return degree + 1;
}

/// is_valid_degree — degree must be non-negative
[[nodiscard]] constexpr bool is_valid_degree(int p) noexcept {
    return p >= 0;
}

// -----------------------------------------------------------------------------
// Min / max helpers for algorithm initialisation
// -----------------------------------------------------------------------------

/// set min to a if a < min, using Tolerance for floating-point safety
template <NumericScalar T>
void set_min(T& min_val, T candidate, const Tolerance<T>& tol) {
    if (tol.lt(candidate, min_val)) min_val = candidate;
}

/// set max to b if b > max, using Tolerance for floating-point safety
template <NumericScalar T>
void set_max(T& max_val, T candidate, const Tolerance<T>& tol) {
    if (tol.gt(candidate, max_val)) max_val = candidate;
}

// -----------------------------------------------------------------------------
// Algorithm result types — these are returned by higher-level algorithms
// and hold intermediate state that multiple algorithms share.
// -----------------------------------------------------------------------------

/**
 * BasisFunctionValues — result of computing p+1 non-zero basis functions
 * at a parameter u.
 */
template <NumericScalar T>
struct BasisFunctionValues {
    std::size_t span;      // knot span index k
    int degree;            // p
    std::vector<T> values;  // N_{k-p,p}(u) … N_{k,p}(u), size = p+1
};

} // namespace nurbs::core

// utilities.hpp — Index helpers, span queries, and algorithmic utilities
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

// types.hpp must come first (defines Point, NURBSPoint, KnotVector, WeightVector).
// types.hpp transitively includes numeric.hpp (NumericScalar_, Tolerance).
#include "types.hpp"

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Span / knot-span queries  (Algorithm A2.1 / A2.2 in NURBS Book)
// -----------------------------------------------------------------------------

/**
 * find_span — Algorithm A2.1 (NURBS Book, 2nd Ed., p.68)
 *
 * Finds the knot span index k such that u ∈ [u_k, u_{k+1}).
 * Special boundary: if u == u_m (last knot), returns n.
 *
 * Performance: O(log n) via binary search.
 */
template <NumericScalar_ T>
std::size_t find_span(std::size_t n, int p, T u, const KnotVector<T>& U) {
    const std::size_t m = U.size() - 1;  // m = n + p + 1

    // Special boundary: at or past the last knot → last span
    if (u >= U[m]) return n;
    // Left boundary: at or before first interior knot span → first span
    if (u <= U[p]) return p;

    std::size_t low  = p;
    std::size_t high = m - p - 1;  // = n

    while (low < high) {
        std::size_t mid = (low + high) / 2;
        if (u >= U[mid]) {
            low = mid + 1;
        } else {
            high = mid;
        }
    }

    return low - 1;
}

/**
 * find_span_linear — O(n) linear scan version for validation / small inputs.
 */
template <NumericScalar_ T>
std::size_t find_span_linear(std::size_t n, int p, T u, const KnotVector<T>& U) {
    const std::size_t m = U.size() - 1;
    if (u >= U[m]) return n;
    for (std::size_t k = p; k < m - p; ++k) {
        if (U[k] <= u && u < U[k + 1]) return k;
    }
    return p;
}

// -----------------------------------------------------------------------------
// Basis function evaluation (Algorithm A2.2 combined with A2.4)
// -----------------------------------------------------------------------------

/**
 * compute_basis_functions — evaluate p+1 non-zero basis functions at u.
 *
 * Fills `b` with N_{k-p,p}(u) … N_{k,p}(u), returns k (the span index).
 */
template <NumericScalar_ T>
std::size_t compute_basis_functions(std::size_t n, int p, T u,
                                    const KnotVector<T>& U,
                                    std::vector<T>& b) {
    std::size_t k = find_span(n, p, u, U);

    b.assign(p + 1, T{0});
    b[0] = T{1};

    if (p == 0) return k;

    std::vector<T> left(p + 1);
    std::vector<T> right(p + 1);
    T saved = T{0};
    T temp  = T{0};

    for (int j = 1; j <= p; ++j) {
        left[j]  = u - U[k + 1 - j];
        right[j] = U[k + j] - u;
        saved    = T{0};

        for (int r = 0; r < j; ++r) {
            T denom = right[r + 1] + left[j - r];
            if (denom == T{0}) {
                b[r] = T{0};
                saved = T{0};
            } else {
                temp      = b[r] / denom;
                b[r]      = saved + right[r + 1] * temp;
                saved     = left[j - r] * temp;
            }
        }
        if (j > 0 && (right[j] + left[0]) == T{0}) {
            b[j] = T{0};
        } else {
            b[j] = saved;
        }
    }

    return k;
}

// -----------------------------------------------------------------------------
// Interval / span utilities
// -----------------------------------------------------------------------------

/** Length of knot span at index k (U[k+1] - U[k]), or 0 for repeated knots. */
template <NumericScalar_ T>
[[nodiscard]] T span_length(std::size_t k, const KnotVector<T>& U) {
    return U[k + 1] - U[k];
}

/** True when U[k] == U[k+1] (span of zero length). */
template <NumericScalar_ T>
[[nodiscard]] bool is_multi_knot(std::size_t k, const KnotVector<T>& U) {
    return U[k + 1] == U[k];
}

/** Number of knot spans with non-zero length. */
template <NumericScalar_ T>
[[nodiscard]] std::size_t num_nonzero_spans(const KnotVector<T>& U) {
    std::size_t count = 0;
    for (std::size_t i = 0; i + 1 < U.size(); ++i)
        if (U[i + 1] != U[i]) ++count;
    return count;
}

// -----------------------------------------------------------------------------
// Domain validation
// -----------------------------------------------------------------------------

/**
 * is_valid_knot_vector — checks structural invariants required for a valid
 * B-spline / NURBS knot vector.
 *
 * Requirements:
 *   1. size() >= 2 * (p + 1)
 *   2. non-decreasing (enforced by KnotVector constructor)
 *   3. for clamped curves: first p+1 knots equal, last p+1 knots equal
 */
template <NumericScalar_ T>
[[nodiscard]] bool is_valid_knot_vector(const KnotVector<T>& U, int p) {
    if (U.size() < static_cast<std::size_t>(2 * (p + 1))) return false;
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
 */
template <NumericScalar_ T>
[[nodiscard]] bool validate_control_point_grid(
    const std::vector<std::vector<NURBSPoint<T>>>& ctrl,
    std::size_t nu,
    std::size_t nv) noexcept
{
    if (ctrl.size() != nv) return false;
    for (std::size_t i = 0; i < nv; ++i)
        if (ctrl[i].size() != nu) return false;
    return true;
}

// -----------------------------------------------------------------------------
// Indexing helpers for 2D surface control point grids
// -----------------------------------------------------------------------------

/** Flat index for surface control point grid at (i, j). */
[[nodiscard]] constexpr std::size_t
surface_control_point_index(std::size_t i, std::size_t j, std::size_t nu) noexcept {
    return i * nu + j;
}

/** Bounds-checked flat index. */
template <NumericScalar_ T>
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

/** polynomial_order = degree + 1 */
[[nodiscard]] constexpr int polynomial_order(int degree) noexcept {
    return degree + 1;
}

/** degree must be non-negative */
[[nodiscard]] constexpr bool is_valid_degree(int p) noexcept {
    return p >= 0;
}

// -----------------------------------------------------------------------------
// Algorithm result types
// -----------------------------------------------------------------------------

/**
 * BasisFunctionValues — result of computing p+1 non-zero basis functions
 * at a parameter u.
 */
template <NumericScalar_ T>
struct BasisFunctionValues {
    std::size_t span;       // knot span index k
    int         degree;     // p
    std::vector<T> values;  // N_{k-p,p}(u) … N_{k,p}(u), size = p+1
};

} // namespace nurbs::core
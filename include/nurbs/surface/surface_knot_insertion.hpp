// surface_knot_insertion.hpp — Surface knot insertion (The NURBS Book, Ch6)
// Algorithm A6.3: insert a knot into a NURBS surface along u or v direction
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/knot_insertion.hpp"

namespace nurbs::surface {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::find_span;

// -----------------------------------------------------------------------------
// Algorithm A6.3 — Surface knot insertion
// -----------------------------------------------------------------------------

/**
 * SurfaceKnotInsertionResult — output of surface_knot_insertion.
 */
template <NumericScalar_ T>
struct SurfaceKnotInsertionResult {
    KnotVector<T> knot_vector_u;                      // U_bar (updated u knot vector)
    KnotVector<T> knot_vector_v;                      // V (unchanged)
    std::vector<std::vector<NURBSPoint<T>>> control_points; // updated control points
};

/**
 * surface_knot_insertion_u — Algorithm A6.3 (NURBS Book, 2nd ed., p.233)
 *
 * Inserts a knot value `u` into a NURBS surface along the u-direction,
 * producing an equivalent surface with one additional control point
 * column.
 *
 * @param u       knot value to insert in u direction
 * @param p_u     polynomial degree in u
 * @param p_v     polynomial degree in v
 * @param U       original knot vector in u
 * @param V       knot vector in v (unchanged)
 * @param Pw      original control points [nv][nu] in homogeneous form
 *
 * @return SurfaceKnotInsertionResult with updated U and control points.
 *
 * The v-direction data (V and control point rows) are unchanged.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceKnotInsertionResult<T>
surface_knot_insertion_u(T u, int p_u, int p_v,
                        const KnotVector<T>& U,
                        const KnotVector<T>& V,
                        const std::vector<std::vector<NURBSPoint<T>>>& Pw) {
    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();
    const std::size_t n = nu - 1;
    const std::size_t m = nv - 1;

    // Find span k for u in U
    std::size_t k = find_span(n, p_u, u, U);

    // Count multiplicity of u in U
    int s = 0;
    for (std::size_t i = 0; i < U.size(); ++i) {
        if (U[i] == u) { s = static_cast<int>(U.multiplicity(i)); break; }
    }

    // At max multiplicity: nothing to do
    if (s == p_u + 1) {
        SurfaceKnotInsertionResult<T> r;
        r.knot_vector_u = U;
        r.knot_vector_v = V;
        r.control_points = Pw;
        return r;
    }

    // Build new u knot vector: insert u at position k
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(U.size() + 1);
    for (std::size_t i = 0; i <= k; ++i) U_bar_vec.push_back(U[i]);
    U_bar_vec.push_back(u);
    for (std::size_t i = k + 1; i < U.size(); ++i) U_bar_vec.push_back(U[i]);
    KnotVector<T> U_bar(std::move(U_bar_vec));

    // New control point grid: nu+1 columns, nv rows
    std::vector<std::vector<NURBSPoint<T>>> Pw_bar(nv, std::vector<NURBSPoint<T>>(nu + 1));

    // Copy unchanged rows for first (k-p_u+1) and last (n-k+1) columns
    // Copy the first (k-p_u+1) columns unchanged
    for (std::size_t j = 0; j <= k - p_u; ++j) {
        for (std::size_t i = 0; i < nv; ++i) {
            Pw_bar[i][j] = Pw[i][j];
        }
    }

    // Copy the last (n - k + 1) columns unchanged
    for (std::size_t j = k + 1; j <= n; ++j) {
        for (std::size_t i = 0; i < nv; ++i) {
            Pw_bar[i][j + 1] = Pw[i][j];
        }
    }

    // Compute new interior control point columns
    int r = p_u - s;
    std::size_t L = k - p_u + 1;
    for (int j = 1; j <= r; ++j) {
        for (std::size_t i = 0; i < nv; ++i) {
            T alpha = (u - U[L + j - 1]) / (U[L + p_u + j - 1] - U[L + j - 1]);
            T alpha_bar = T{1} - alpha;
            Pw_bar[i][L + j - 1] = Pw[i][L + j - 1] * alpha_bar + Pw_bar[i][L + j - 2] * alpha;
        }
    }

    SurfaceKnotInsertionResult<T> result;
    result.knot_vector_u = std::move(U_bar);
    result.knot_vector_v = V;
    result.control_points = std::move(Pw_bar);
    return result;
}

/**
 * surface_knot_insertion_v — Algorithm A6.3 (NURBS Book, 2nd ed., p.233)
 *
 * Inserts a knot value `v` into a NURBS surface along the v-direction,
 * producing an equivalent surface with one additional control point row.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceKnotInsertionResult<T>
surface_knot_insertion_v(T v, int p_u, int p_v,
                        const KnotVector<T>& U,
                        const KnotVector<T>& V,
                        const std::vector<std::vector<NURBSPoint<T>>>& Pw) {
    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();
    const std::size_t n = nu - 1;
    const std::size_t m = nv - 1;

    // Find span l for v in V
    std::size_t l = find_span(m, p_v, v, V);

    // Count multiplicity of v in V
    int s = 0;
    for (std::size_t i = 0; i < V.size(); ++i) {
        if (V[i] == v) { s = static_cast<int>(V.multiplicity(i)); break; }
    }

    if (s == p_v + 1) {
        SurfaceKnotInsertionResult<T> r;
        r.knot_vector_u = U;
        r.knot_vector_v = V;
        r.control_points = Pw;
        return r;
    }

    // Build new v knot vector
    std::vector<T> V_bar_vec;
    V_bar_vec.reserve(V.size() + 1);
    for (std::size_t i = 0; i <= l; ++i) V_bar_vec.push_back(V[i]);
    V_bar_vec.push_back(v);
    for (std::size_t i = l + 1; i < V.size(); ++i) V_bar_vec.push_back(V[i]);
    KnotVector<T> V_bar(std::move(V_bar_vec));

    // New control point grid: nv+1 rows, nu columns
    std::vector<std::vector<NURBSPoint<T>>> Pw_bar(nv + 1, std::vector<NURBSPoint<T>>(nu));

    // Copy the first (l-p_v+1) rows unchanged
    for (std::size_t i = 0; i <= l - p_v; ++i) {
        for (std::size_t j = 0; j < nu; ++j) {
            Pw_bar[i][j] = Pw[i][j];
        }
    }

    // Copy the last (m - l + 1) rows unchanged
    for (std::size_t i = l + 1; i <= m; ++i) {
        for (std::size_t j = 0; j < nu; ++j) {
            Pw_bar[i + 1][j] = Pw[i][j];
        }
    }

    // Compute new interior control point rows
    int r = p_v - s;
    std::size_t L = l - p_v + 1;
    for (int j = 1; j <= r; ++j) {
        for (std::size_t jj = 0; jj < nu; ++jj) {
            T alpha = (v - V[L + j - 1]) / (V[L + p_v + j - 1] - V[L + j - 1]);
            T alpha_bar = T{1} - alpha;
            Pw_bar[L + j - 1][jj] = Pw[L + j - 1][jj] * alpha_bar + Pw_bar[L + j - 2][jj] * alpha;
        }
    }

    SurfaceKnotInsertionResult<T> result;
    result.knot_vector_u = U;
    result.knot_vector_v = std::move(V_bar);
    result.control_points = std::move(Pw_bar);
    return result;
}

/**
 * surface_knot_insertion — general overload that inserts in u by default.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceKnotInsertionResult<T>
surface_knot_insertion(T u, int p_u, int p_v,
                      const KnotVector<T>& U,
                      const KnotVector<T>& V,
                      const std::vector<std::vector<NURBSPoint<T>>>& Pw,
                      bool direction_u = true) {
    if (direction_u)
        return surface_knot_insertion_u(u, p_u, p_v, U, V, Pw);
    else
        return surface_knot_insertion_v(u, p_u, p_v, U, V, Pw);
}

} // namespace nurbs::surface
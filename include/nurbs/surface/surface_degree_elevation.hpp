// surface_degree_elevation.hpp — Surface degree elevation (The NURBS Book, Ch6)
// Algorithm A6.4: elevate the degree of a NURBS surface in u or v direction
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/degree_elevation.hpp"

namespace nurbs::surface {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;

// -----------------------------------------------------------------------------
// Algorithm A6.4 — Surface degree elevation
// -----------------------------------------------------------------------------

/**
 * SurfaceDegreeElevationResult — output of surface_degree_elevation.
 */
template <NumericScalar_ T>
struct SurfaceDegreeElevationResult {
    int degree_u;  // p_bar (= p + t_u)
    int degree_v;  // q_bar (= q + t_v)
    KnotVector<T> knot_vector_u;   // U_bar
    KnotVector<T> knot_vector_v;   // V_bar
    std::vector<std::vector<NURBSPoint<T>>> control_points; // updated control points
};

/**
 * surface_degree_elevation_u — Algorithm A6.4 (NURBS Book, 2nd ed., p.236)
 *
 * Elevates a NURBS surface by t degrees in the u-direction only.
 * The v-direction is unchanged.
 *
 * @param p_u    original polynomial degree in u
 * @param p_v    polynomial degree in v (unchanged)
 * @param U      original knot vector in u
 * @param V      knot vector in v (unchanged)
 * @param Pw     original control points [nv][nu] in homogeneous form
 * @param t      number of degrees to elevate in u (t >= 1)
 *
 * @return SurfaceDegreeElevationResult with elevated u-degree surface data.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceDegreeElevationResult<T>
surface_degree_elevation_u(int p_u, int p_v,
                         const KnotVector<T>& U,
                         const KnotVector<T>& V,
                         const std::vector<std::vector<NURBSPoint<T>>>& Pw,
                         int t) {
    if (t < 1) throw std::invalid_argument("surface_degree_elevation_u: t must be >= 1");
    if (p_u < 0) throw std::invalid_argument("surface_degree_elevation_u: degree must be >= 0");

    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();

    // Elevate degree in u for each v-row of control points
    // New nu after elevation: nu + t
    const std::size_t nu_bar = nu + static_cast<std::size_t>(t);
    std::vector<std::vector<NURBSPoint<T>>> Pw_bar(nv, std::vector<NURBSPoint<T>>(nu_bar));

    // Use basis-level degree elevation (Algorithm A3.4) on each v=row
    for (std::size_t i = 0; i < nv; ++i) {
        // Extract homogeneous control points for row i
        std::vector<NURBSPoint<T>> row = Pw[i];
        auto row_result = nurbs::basis::elevate_degree(U, p_u, row, t);

        // Copy elevated row control points
        for (std::size_t j = 0; j < row_result.control_points.size(); ++j)
            Pw_bar[i][j] = row_result.control_points[j];
    }

    // Build the elevated knot vector U_bar:
    // - Copy first (p_u+1) and last (p_u+1) knots unchanged
    // - Insert t copies of each interior knot
    std::vector<T> U_bar_vec;
    U_bar_vec.reserve(U.size() + static_cast<std::size_t>(t) * U.num_interior());

    for (int i = 0; i <= p_u; ++i) U_bar_vec.push_back(U[static_cast<std::size_t>(i)]);
    for (std::size_t i = static_cast<std::size_t>(p_u + 1);
         i + p_u + 1 < U.size() - 1; ++i) {
        U_bar_vec.push_back(U[i]);
        if (U[i] != U[i + 1]) { // interior knot
            for (int j = 0; j < t; ++j) U_bar_vec.push_back(U[i]);
        }
    }
    for (int i = static_cast<int>(U.size()) - p_u - 1; i < static_cast<int>(U.size()); ++i)
        U_bar_vec.push_back(U[static_cast<std::size_t>(i)]);

    SurfaceDegreeElevationResult<T> result;
    result.degree_u = p_u + t;
    result.degree_v = p_v;
    result.knot_vector_u = KnotVector<T>(std::move(U_bar_vec));
    result.knot_vector_v = V;
    result.control_points = std::move(Pw_bar);
    return result;
}

/**
 * surface_degree_elevation_v — Algorithm A6.4 (NURBS Book, 2nd ed., p.236)
 *
 * Elevates a NURBS surface by t degrees in the v-direction only.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceDegreeElevationResult<T>
surface_degree_elevation_v(int p_u, int p_v,
                         const KnotVector<T>& U,
                         const KnotVector<T>& V,
                         const std::vector<std::vector<NURBSPoint<T>>>& Pw,
                         int t) {
    if (t < 1) throw std::invalid_argument("surface_degree_elevation_v: t must be >= 1");
    if (p_v < 0) throw std::invalid_argument("surface_degree_elevation_v: degree must be >= 0");

    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();
    const std::size_t nv_bar = nv + static_cast<std::size_t>(t);

    std::vector<std::vector<NURBSPoint<T>>> Pw_bar(nv_bar, std::vector<NURBSPoint<T>>(nu));

    // Elevate degree in v for each u-column of control points
    for (std::size_t j = 0; j < nu; ++j) {
        std::vector<NURBSPoint<T>> col(nv);
        for (std::size_t i = 0; i < nv; ++i)
            col[i] = Pw[i][j];

        auto col_result = nurbs::basis::elevate_degree(V, p_v, col, t);

        for (std::size_t i = 0; i < col_result.control_points.size(); ++i)
            Pw_bar[i][j] = col_result.control_points[i];
    }

    // Build elevated V knot vector
    std::vector<T> V_bar_vec;
    V_bar_vec.reserve(V.size() + static_cast<std::size_t>(t) * V.num_interior());

    for (int i = 0; i <= p_v; ++i) V_bar_vec.push_back(V[static_cast<std::size_t>(i)]);
    for (std::size_t i = static_cast<std::size_t>(p_v + 1);
         i + p_v + 1 < V.size() - 1; ++i) {
        V_bar_vec.push_back(V[i]);
        if (V[i] != V[i + 1]) {
            for (int j = 0; j < t; ++j) V_bar_vec.push_back(V[i]);
        }
    }
    for (int i = static_cast<int>(V.size()) - p_v - 1; i < static_cast<int>(V.size()); ++i)
        V_bar_vec.push_back(V[static_cast<std::size_t>(i)]);

    SurfaceDegreeElevationResult<T> result;
    result.degree_u = p_u;
    result.degree_v = p_v + t;
    result.knot_vector_u = U;
    result.knot_vector_v = KnotVector<T>(std::move(V_bar_vec));
    result.control_points = std::move(Pw_bar);
    return result;
}

/**
 * surface_degree_elevation — elevate in both u and v directions.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceDegreeElevationResult<T>
surface_degree_elevation(int p_u, int p_v,
                       const KnotVector<T>& U,
                       const KnotVector<T>& V,
                       const std::vector<std::vector<NURBSPoint<T>>>& Pw,
                       int t_u, int t_v) {
    // First elevate in u, then in v
    auto u_result = surface_degree_elevation_u(p_u, p_v, U, V, Pw, t_u);
    auto v_result = surface_degree_elevation_v(
        u_result.degree_u, u_result.degree_v,
        u_result.knot_vector_u, u_result.knot_vector_v,
        u_result.control_points, t_v);
    return v_result;
}

} // namespace nurbs::surface
// construct_surface.hpp — NURBS surface construction (The NURBS Book, Ch6)
// Algorithm A6.1: construct a NURBS surface from control points and knot vectors
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/bspline_basis.hpp"
#include "../basis/knot_insertion.hpp"
#include "nurbs_surface.hpp"

namespace nurbs::surface {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::WeightVector;
using nurbs::core::find_span;
using nurbs::core::compute_basis_functions;
using nurbs::core::validate_control_point_grid;

// -----------------------------------------------------------------------------
// Algorithm A6.1 — Surface construction
// -----------------------------------------------------------------------------

/**
 * SurfaceConstructionResult — output of construct_surface.
 */
template <NumericScalar_ T>
struct SurfaceConstructionResult {
    int degree_u;
    int degree_v;
    KnotVector<T> knot_vector_u;
    KnotVector<T> knot_vector_v;
    std::vector<std::vector<NURBSPoint<T>>> control_points; // size [nv][nu]
};

/**
 * construct_surface — Algorithm A6.1 (NURBS Book, 2nd ed., p.228)
 *
 * Validates and wraps a NURBS surface from its fundamental data.
 * A NURBS surface of bi-degree (p, q) with (n+1)×(m+1) control points
 * requires knot vectors U (size n+p+2) and V (size m+q+2).
 *
 * @param p_u     polynomial degree in u direction
 * @param p_v     polynomial degree in v direction
 * @param U       knot vector in u direction
 * @param V       knot vector in v direction
 * @param Pw      2D grid of control points in homogeneous form [m+1][n+1]
 *
 * @return SurfaceConstructionResult wrapping the surface data.
 *
 * Algorithm A6.1 mainly provides structural validation.  Evaluation uses
 * the homogeneous tensor-product De Boor algorithm.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceConstructionResult<T>
construct_surface(int p_u, int p_v,
                 const KnotVector<T>& U,
                 const KnotVector<T>& V,
                 const std::vector<std::vector<NURBSPoint<T>>>& Pw) {
    if (p_u < 0) throw std::invalid_argument("construct_surface: p_u must be >= 0");
    if (p_v < 0) throw std::invalid_argument("construct_surface: p_v must be >= 0");

    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();

    // Validate: Pw must be rectangular
    if (!Pw.empty() && !validate_control_point_grid(Pw, nu, nv))
        throw std::invalid_argument("construct_surface: Pw must be rectangular [nv][nu]");

    // n+1 control points in u direction, m+1 in v direction
    const std::size_t n = nu - 1;
    const std::size_t m = nv - 1;

    // Knot vector size invariants: m_u = n + p_u + 1, m_v = m + q + 1
    if (U.size() != n + static_cast<std::size_t>(p_u) + 2)
        throw std::invalid_argument("construct_surface: U size must equal n + p_u + 2");
    if (V.size() != m + static_cast<std::size_t>(p_v) + 2)
        throw std::invalid_argument("construct_surface: V size must equal m + p_v + 2");

    SurfaceConstructionResult<T> result;
    result.degree_u = p_u;
    result.degree_v = p_v;
    result.knot_vector_u = U;
    result.knot_vector_v = V;
    result.control_points = Pw;
    return result;
}

/**
 * construct_surface_from_cartesian — build a NURBS surface from Cartesian points.
 *
 * @param p_u    polynomial degree in u
 * @param p_v    polynomial degree in v
 * @param U      knot vector in u
 * @param V      knot vector in v
 * @param P      Cartesian control points (flat, row-major or column-major)
 * @param nu     number of control points in u direction
 * @param nv     number of control points in v direction
 * @param W      weights (uniform if omitted)
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceConstructionResult<T>
construct_surface_from_cartesian(int p_u, int p_v,
                                const KnotVector<T>& U,
                                const KnotVector<T>& V,
                                const std::vector<nurbs::core::Point<T>>& P,
                                std::size_t nu, std::size_t nv,
                                const WeightVector<T>& W) {
    if (P.size() != nu * nv)
        throw std::invalid_argument("construct_surface_from_cartesian: P size must equal nu*nv");

    std::vector<std::vector<NURBSPoint<T>>> Pw(nv, std::vector<NURBSPoint<T>>(nu));
    for (std::size_t i = 0; i < nv; ++i) {
        for (std::size_t j = 0; j < nu; ++j) {
            std::size_t idx = i * nu + j;
            T w = W.empty() ? T{1} : W[idx];
            Pw[i][j] = NURBSPoint<T>(P[idx], w);
        }
    }

    return construct_surface(p_u, p_v, U, V, Pw);
}

} // namespace nurbs::surface
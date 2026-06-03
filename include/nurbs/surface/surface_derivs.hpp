// surface_derivs.hpp — NURBS surface derivative evaluation (The NURBS Book, Ch6)
// Algorithm A6.2: evaluate a NURBS surface point and its parametric partial derivatives
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/bspline_basis.hpp"
#include "../curve/curve_derivatives.hpp"

namespace nurbs::surface {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::find_span;
using nurbs::core::compute_basis_functions;

// -----------------------------------------------------------------------------
// Algorithm A6.2 — Surface derivatives
// -----------------------------------------------------------------------------

/**
 * SurfaceDerivatives — result bundle from surface derivative evaluation.
 */
template <NumericScalar_ T>
struct SurfaceDerivatives {
    int degree_u; // p
    int degree_v; // q
    //derivs[d_u][d_v] is a 2D array of NURBSPoint for mixed partial derivative
    //of order (d_u, d_v). Storage: flattened row-major.
    std::vector<std::vector<NURBSPoint<T>>> derivs;
};

/**
 * surface_derivatives — Algorithm A6.2 (NURBS Book, 2nd ed., p.230)
 *
 * Evaluates a NURBS surface and its parametric partial derivatives up to
 * orders d_u in u and d_v in v at a given parameter point (u, v).
 *
 * The algorithm uses the tensor-product structure of NURBS surfaces:
 *   S(u,v) = Σ_{i=0}^{n} Σ_{j=0}^{m} R_{i,j}(u,v) · P_{i,j}
 * where the rational basis functions R_{i,j} = N_{i,p}(u)·N_{j,q}(v)·w_{i,j} / w(u,v)
 *
 * @param u       parameter value in u direction
 * @param v       parameter value in v direction
 * @param p       polynomial degree in u
 * @param q       polynomial degree in v
 * @param U       knot vector in u
 * @param V       knot vector in v
 * @param Pw      control points in homogeneous form (size nv × nu)
 * @param d_u     number of u-derivatives to compute (0 = value only)
 * @param d_v     number of v-derivatives to compute (0 = value only)
 *
 * @return SurfaceDerivatives with a [d_u+1] × [d_v+1] array of NURBSPoint
 *         representing S(u,v) and its partial derivatives.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceDerivatives<T>
surface_derivatives(T u, T v,
                   int p, int q,
                   const KnotVector<T>& U,
                   const KnotVector<T>& V,
                   const std::vector<std::vector<NURBSPoint<T>>>& Pw,
                   int d_u, int d_v) {
    const std::size_t nu = Pw.empty() ? 0 : Pw[0].size();
    const std::size_t nv = Pw.size();
    const std::size_t n = nu - 1;
    const std::size_t m = nv - 1;

    const int d_u_min = std::min(d_u, p);
    const int d_v_min = std::min(d_v, q);

    SurfaceDerivatives<T> result;
    result.degree_u = p;
    result.degree_v = q;
    result.derivs.resize(static_cast<std::size_t>(d_u_min) + 1);
    for (int du = 0; du <= d_u_min; ++du) {
        result.derivs[static_cast<std::size_t>(du)].resize(
            static_cast<std::size_t>(d_v_min) + 1);
    }

    // ---- Find spans k and l ----
    std::size_t k = find_span(n, p, u, U);
    std::size_t l = find_span(m, q, v, V);

    // ---- Compute basis functions in u and v ----
    std::vector<T> b_u(p + 1);
    std::vector<T> b_v(q + 1);
    compute_basis_functions(n, p, u, U, b_u);
    compute_basis_functions(m, q, v, V, b_v);

    // ---- Temporary arrays for the tensor-product algorithm ----
    // We use the homogeneous De Boor approach similar to curves, but
    // tensor-product: first evaluate the v-direction control points,
    // then use them as input for u-direction De Boor.
    //
    // Algorithm A6.2 (simplified):
    // 1. For each v-level j (0..q), compute curve derivatives in u using
    //    the control point column Pw[*, j] to get CK_u[du][i]
    // 2. Use these as control points for v-direction derivative evaluation.

    // Temporary: store homogeneous control points for each v-basis column
    // Qw_temp[i][j] = control point at (i,j) in homogeneous coordinates
    std::vector<std::vector<NURBSPoint<T>>> Qw_temp(nu, std::vector<NURBSPoint<T>>(nv));
    for (std::size_t i = 0; i < nu; ++i)
        for (std::size_t j = 0; j < nv; ++j)
            Qw_temp[i][j] = Pw[j][i]; // Note: Pw is [nv][nu], transposed here

    // Step 1: compute u-derivatives for each of the q+1 v-columns
    // For each l-v-span (v-direction basis function index), build a curve
    // and compute derivatives.
    for (int w = 0; w <= q; ++w) {
        // Extract the w-th column of control points for the v basis function
        std::vector<NURBSPoint<T>> col(nu);
        for (std::size_t i = 0; i < nu; ++i)
            col[i] = Pw[l - q + w][i]; // Pw is [nv][nu]

        // Compute curve derivatives in u at parameter u for this column
        // using Algorithm A5.2 (curves)
        auto cders = nurbs::curve::curve_derivatives(u, p, U, col, d_u_min);

        // Store the du-th derivative into the result array
        for (int du = 0; du <= d_u_min; ++du) {
            result.derivs[static_cast<std::size_t>(du)]
                         [static_cast<std::size_t>(w)] = cders.values[static_cast<std::size_t>(du)];
        }
    }

    // Step 2: now do v-direction derivative evaluation using the
    // results from step 1 as control points
    // For each u-derivative order du, treat it as a curve along v and
    // compute v-derivatives.

    std::vector<NURBSPoint<T>> temp_v(q + 1);
    for (int du = 0; du <= d_u_min; ++du) {
        // Use the u-derivative results as v-direction "control points"
        // Extract column from result.derivs[du][w]
        std::vector<NURBSPoint<T>> v_ctrl(q + 1);
        for (int w = 0; w <= q; ++w) {
            v_ctrl[static_cast<std::size_t>(w)] =
                result.derivs[static_cast<std::size_t>(du)]
                             [static_cast<std::size_t>(w)];
        }

        // Evaluate v-direction derivatives using Algorithm A5.2
        auto vders = nurbs::curve::curve_derivatives(v, q, V, v_ctrl, d_v_min);

        for (int dv = 0; dv <= d_v_min; ++dv) {
            result.derivs[static_cast<std::size_t>(du)]
                         [static_cast<std::size_t>(dv)] = vders.values[static_cast<std::size_t>(dv)];
        }
    }

    return result;
}

/**
 * surface_derivatives — convenience overload with default derivative order = degree.
 */
template <NumericScalar_ T>
[[nodiscard]] SurfaceDerivatives<T>
surface_derivatives(T u, T v,
                   int p, int q,
                   const KnotVector<T>& U,
                   const KnotVector<T>& V,
                   const std::vector<std::vector<NURBSPoint<T>>>& Pw) {
    return surface_derivatives(u, v, p, q, U, V, Pw, p, q);
}

} // namespace nurbs::surface
// curve_derivatives.hpp — NURBS curve derivative evaluation (The NURBS Book, Ch5)
// Algorithm A5.2: evaluate curve and its parametric derivatives
#pragma once

#include <algorithm>
#include <cstddef>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/bspline_basis.hpp"
#include "nurbs_curve.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::find_span;
using nurbs::core::compute_basis_functions;

// -----------------------------------------------------------------------------
// Algorithm A5.2 — Curve derivatives
// -----------------------------------------------------------------------------

/**
 * CurveDerivatives — result bundle from derivative evaluation.
 */
template <NumericScalar_ T>
struct CurveDerivatives {
    int degree;                                  // curve degree p
    std::vector<NURBSPoint<T>> values;           // CK[0..d]: point + derivatives in homogeneous form
};

/**
 * curve_derivatives — Algorithm A5.2 (NURBS Book, 2nd ed., p.213)
 *
 * Computes the point on the curve and its parametric derivatives up to order d
 * at a given parameter value u.
 *
 * @param u    parameter value
 * @param p    polynomial degree
 * @param U    knot vector
 * @param Pw   control points in homogeneous coordinates
 * @param d    number of derivatives to compute (0 = point only, 1 = point+1st deriv, …)
 *
 * @return CurveDerivatives with (d+1) homogeneous points: C(u), C'(u), …, C^{(d)}(u)
 *
 * Algorithm A5.2 works in homogeneous coordinates and then dehomogenizes.
 */
template <NumericScalar_ T>
[[nodiscard]] CurveDerivatives<T>
curve_derivatives(T u, int p, const KnotVector<T>& U,
                  const std::vector<NURBSPoint<T>>& Pw,
                  int d) {
    const std::size_t n = Pw.size() - 1;
    const std::size_t m = U.size() - 1;
    const int dmin = std::min(d, p);

    CurveDerivatives<T> result;
    result.degree = p;
    result.values.resize(dmin + 1);

    if (d < 0) return result;

    // Find span
    std::size_t k = find_span(n, p, u, U);

    // ---- Compute all basis function derivatives N_{i,p}^{(r)}(u) for r=0..d, i=0..p ----
    std::vector<std::vector<T>> ders(dmin + 1, std::vector<T>(p + 1, T{0}));
    {
        std::vector<std::vector<T>> ndu(p + 1, std::vector<T>(p + 1, T{0}));
        std::vector<T> left(p + 1, T{0});
        std::vector<T> right(p + 1, T{0});

        ndu[0][0] = T{1};
        for (int j = 1; j <= p; ++j) {
            left[j]  = u - U[k + 1 - j];
            right[j] = U[k + j] - u;
            T saved  = T{0};
            for (int r = 0; r < j; ++r) {
                ndu[j][r] = right[r + 1] + left[j - r];
                T temp = (ndu[j][r] == T{0})
                         ? T{0}
                         : (saved + right[r + 1] * ders[0][r] / ndu[j][r]);
                ders[0][r] = temp;
                saved = left[j - r] * temp;
            }
            ders[0][j] = saved;
        }

        // Derivative table computation (Algorithm A2.3 simplified)
        for (int r = 1; r <= dmin; ++r) {
            for (int i = 0; i <= p; ++i) {
                T sum = T{0};
                int ri = r - 1;
                if (p - ri >= 0) {
                    sum = (ders[ri][i] * static_cast<T>(r)
                           - ders[ri][i - 1]) / (U[k + p + ri - i + 1] - U[k + ri - i + 1]);
                }
                ders[r][i] = sum;
            }
        }
    }

    // ---- Compute curve derivatives CK[r] = Σ_{i=0}^{n} N_{i,p}^{(r)}(u) · Pw_i ----
    for (int r = 0; r <= dmin; ++r) {
        NURBSPoint<T> sum;
        for (int i = 0; i <= p; ++i) {
            std::size_t idx = k - p + static_cast<std::size_t>(i);
            if (idx < Pw.size()) {
                sum = sum + Pw[idx] * ders[r][i];
            }
        }
        result.values[static_cast<std::size_t>(r)] = sum;
    }

    // Dehomogenize if needed for display purposes (return homogeneous for exact arithmetic)
    return result;
}

/**
 * curve_derivatives — convenience wrapper that accepts a NURBSCurve.
 */
template <NumericScalar_ T>
[[nodiscard]] CurveDerivatives<T>
curve_derivatives(T u, int p, const KnotVector<T>& U,
                  const std::vector<NURBSPoint<T>>& Pw) {
    return curve_derivatives(u, p, U, Pw, p);
}

} // namespace nurbs::curve
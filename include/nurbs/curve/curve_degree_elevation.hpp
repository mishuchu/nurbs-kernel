// curve_degree_elevation.hpp — Curve degree elevation (The NURBS Book, Ch5)
// Algorithm A5.3: elevate the degree of a NURBS curve by t
#pragma once

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/degree_elevation.hpp"
#include "../basis/knot_insertion.hpp"
#include "nurbs_curve.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::WeightVector;

// -----------------------------------------------------------------------------
// Algorithm A5.3 — Curve degree elevation
// -----------------------------------------------------------------------------

/**
 * CurveDegreeElevationResult — output of curve_degree_elevation.
 */
template <NumericScalar_ T>
struct CurveDegreeElevationResult {
    int degree;                                // p_bar (= p + t)
    KnotVector<T> knot_vector;                 // U_bar
    std::vector<NURBSPoint<T>> control_points; // Qw_bar
};

/**
 * curve_degree_elevation — Algorithm A5.3 (NURBS Book, 2nd ed., p.205)
 *
 * Elevates a NURBS curve of degree p to degree p+t without changing its shape.
 * The result is an equivalent curve with more control points and a revised
 * knot vector.
 *
 * This is a wrapper around the basis-level Algorithm A3.4 that handles the
 * NURBS (rational) case by converting to homogeneous form.
 *
 * @param p     original polynomial degree
 * @param U     original knot vector
 * @param Pw    original control points in homogeneous form
 * @param t     number of degrees to elevate (t >= 1)
 *
 * @return CurveDegreeElevationResult with elevated curve data.
 */
template <NumericScalar_ T>
[[nodiscard]] CurveDegreeElevationResult<T>
curve_degree_elevation(int p, const KnotVector<T>& U,
                      const std::vector<NURBSPoint<T>>& Pw,
                      int t) {
    if (t < 1) throw std::invalid_argument("curve_degree_elevation: t must be >= 1");
    if (p < 0) throw std::invalid_argument("curve_degree_elevation: degree must be >= 0");

    // Use the basis-level degree elevation (Algorithm A3.4)
    auto basis_result = nurbs::basis::elevate_degree(U, p, Pw, t);

    CurveDegreeElevationResult<T> result;
    result.degree = p + t;
    result.knot_vector = std::move(basis_result.knot_vector);
    result.control_points = std::move(basis_result.control_points);
    return result;
}

/**
 * curve_elevate_degree — convenience overload that returns a new NURBSCurve.
 *
 * @param p   original degree
 * @param U   original knot vector
 * @param Pw  original homogeneous control points
 * @param t   elevation amount
 *
 * @return NURBSCurve of degree p+t
 */
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_elevate_degree(int p, const KnotVector<T>& U,
                     const std::vector<NURBSPoint<T>>& Pw,
                     int t) {
    auto res = curve_degree_elevation(p, U, Pw, t);
    return NURBSCurve<T>(res.degree, std::move(res.knot_vector),
                        std::move(res.control_points));
}

} // namespace nurbs::curve
// curve_knot_insertion.hpp — Curve knot insertion (The NURBS Book, Ch5)
// Algorithm A5.4: insert a knot into a NURBS curve (wrapper around Algorithm A3.3)
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/knot_insertion.hpp"
#include "nurbs_curve.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;

// -----------------------------------------------------------------------------
// Algorithm A5.4 — Curve knot insertion
// -----------------------------------------------------------------------------

/**
 * CurveKnotInsertionResult — output of curve_knot_insertion.
 */
template <NumericScalar_ T>
struct CurveKnotInsertionResult {
    int degree;                                 // p (unchanged)
    KnotVector<T> knot_vector;                 // U_bar
    std::vector<NURBSPoint<T>> control_points; // Qw_bar
};

/**
 * curve_knot_insertion — Algorithm A5.4 (NURBS Book, 2nd ed., p.151)
 *
 * Inserts a new knot value `u` into a NURBS curve defined by (p, U, Pw),
 * producing an equivalent curve with potentially one additional control point.
 *
 * This is a thin wrapper around the basis-level Algorithm A3.3
 * (knot_insertion.hpp) that preserves the curve API.
 *
 * @param u    knot value to insert (must lie in the parameter domain)
 * @param p    polynomial degree
 * @param U    original knot vector
 * @param Pw   original control points in homogeneous form
 *
 * @return CurveKnotInsertionResult with updated curve data.
 */
template <NumericScalar_ T>
[[nodiscard]] CurveKnotInsertionResult<T>
curve_knot_insertion(T u, int p, const KnotVector<T>& U,
                    const std::vector<NURBSPoint<T>>& Pw) {
    auto res = nurbs::basis::insert_knot(u, p, U, Pw);

    CurveKnotInsertionResult<T> result;
    result.degree = p;
    result.knot_vector = std::move(res.knot_vector);
    result.control_points = std::move(res.control_points);
    return result;
}

/**
 * curve_insert_knot — convenience wrapper that returns a new NURBSCurve.
 *
 * @param u   knot value to insert
 * @param p   polynomial degree
 * @param U   knot vector
 * @param Pw  homogeneous control points
 *
 * @return new NURBSCurve with the knot inserted.
 */
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
curve_insert_knot(T u, int p, const KnotVector<T>& U,
                 const std::vector<NURBSPoint<T>>& Pw) {
    auto res = curve_knot_insertion(u, p, U, Pw);
    return NURBSCurve<T>(res.degree, std::move(res.knot_vector),
                         std::move(res.control_points));
}

} // namespace nurbs::curve
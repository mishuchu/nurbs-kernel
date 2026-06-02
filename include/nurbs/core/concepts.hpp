// concepts.hpp — C++20 concept constraints for NURBS types
#pragma once

#include <concepts>
#include <utility>
#include <vector>

// types.hpp transitively includes numeric.hpp (NumericScalar_, Tolerance, etc.)
// All classes (Point, NURBSPoint, KnotVector, WeightVector) are fully defined
// by the time this header is parsed — no forward declarations needed.
#include "types.hpp"

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Parametric entity concept
// -----------------------------------------------------------------------------

/**
 * ParametricEntity — anything that exposes a valid parametric range.
 */
template <typename E, typename T = double>
concept ParametricEntity = requires(const E& e) {
    { e.parameter_domain() } -> std::same_as<std::pair<T, T>>;
    requires NumericScalar_<T>;
};

// -----------------------------------------------------------------------------
// Basis function support
// -----------------------------------------------------------------------------

/**
 * BasisFunctionProvider — objects that carry a knot vector and polynomial
 * degree and can supply B-spline basis function values.
 */
template <typename B, typename T = double>
concept BasisFunctionProvider = requires(const B& b) {
    { b.degree() } -> std::same_as<int>;
    { b.knot_vector() } -> std::same_as<const KnotVector<T>&>;
    requires NumericScalar_<T>;
    requires b.degree() >= 0;
};

// -----------------------------------------------------------------------------
// Curve concept
// -----------------------------------------------------------------------------

/**
 * CurveConcept — a parametric curve in ℝ^n.
 */
template <typename C, typename T = double>
concept CurveConcept = ParametricEntity<C, T> && BasisFunctionProvider<C, T>
    && requires(const C& c, T u) {
        { c.evaluate(u) } -> std::same_as<NURBSPoint<T>>;
        { c.control_points() } -> std::same_as<std::vector<NURBSPoint<T>>>;
        { c.weights() } -> std::same_as<WeightVector<T>>;
        requires NumericScalar_<T>;
};

// -----------------------------------------------------------------------------
// Surface concept
// -----------------------------------------------------------------------------

/**
 * SurfaceConcept — a bivariate parametric surface in ℝ^n.
 */
template <typename S, typename T = double>
concept SurfaceConcept = requires(const S& s, T u, T v) {
    { s.evaluate(u, v) } -> std::same_as<NURBSPoint<T>>;
    { s.u_degree() } -> std::same_as<int>;
    { s.v_degree() } -> std::same_as<int>;
    { s.u_knot_vector() } -> std::same_as<const KnotVector<T>&>;
    { s.v_knot_vector() } -> std::same_as<const KnotVector<T>&>;
    { s.control_points() } -> std::same_as<std::vector<std::vector<NURBSPoint<T>>>>;
    { s.weights() } -> std::same_as<WeightVector<T>>;
    { s.parameter_domain_u() } -> std::same_as<std::pair<T, T>>;
    { s.parameter_domain_v() } -> std::same_as<std::pair<T, T>>;
    requires NumericScalar_<T>;
    requires s.u_degree() >= 0;
    requires s.v_degree() >= 0;
};

// -----------------------------------------------------------------------------
// Point / weight helpers
// -----------------------------------------------------------------------------

/// Concept for a geometric point or vector in ℝ^n
template <typename P, typename T = double>
concept PointLike = requires(const P& a, const P& b, T s) {
    { a + b } -> std::same_as<P>;
    { a - b } -> std::same_as<P>;
    { a * s } -> std::same_as<P>;
    { s * a } -> std::same_as<P>;
    { a / s } -> std::same_as<P>;
    requires NumericScalar_<T>;
};

/// WeightVector — a range of positive scalars associated with control points
template <typename W, typename T = double>
concept WeightRange = requires(const W& w) {
    { w.size() } -> std::same_as<std::size_t>;
    { w[std::declval<std::size_t>()] } -> std::same_as<T>;
    requires NumericScalar_<T>;
};

} // namespace nurbs::core

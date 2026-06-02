// concepts.hpp — C++20 concept constraints for NURBS types
#pragma once

#include <concepts>
#include <type_traits>

namespace nurbs::core {

// Forward declarations
template <std::floating_point T>
class Point;

template <std::floating_point T>
class NURBSPoint;

template <std::floating_point T>
class KnotVector;

template <std::floating_point T>
class WeightVector;

// -----------------------------------------------------------------------------
// Scalar numeric types
// -----------------------------------------------------------------------------

/// True for builtin floating-point types usable as NURBS coordinate scalars.
template <typename T>
concept Numeric = std::floating_point<T>;

// -----------------------------------------------------------------------------
// Parametric entity concept — anything that exposes a valid parametric range
// -----------------------------------------------------------------------------

/**
 * ParametricEntity
 *
 * Models a parametric geometry object that has a parameter domain.
 * Concrete examples: CurveConcept, SurfaceConcept.
 *
 * Requirements:
 *   - `parameter_domain()` → std::pair<T, T> with domain.first ≤ domain.second
 */
template <typename E, typename T = double>
concept ParametricEntity = requires(const E& e) {
    { e.parameter_domain() } -> std::same_as<std::pair<T, T>>;
    requires Numeric<T>;
};

// -----------------------------------------------------------------------------
// Basis function support — objects that can evaluate B-spline basis functions
// -----------------------------------------------------------------------------

/**
 * BasisFunctionProvider
 *
 * Objects that carry a knot vector and polynomial degree and can therefore
 * supply B-spline basis function values on demand.
 *
 * Requirements:
 *   - `degree()`   → int (≥ 0)
 *   - `knot_vector()` → const KnotVector<T>&
 */
template <typename B, typename T = double>
concept BasisFunctionProvider = requires(const B& b) {
    { b.degree() } -> std::same_as<int>;
    { b.knot_vector() } -> std::same_as<const KnotVector<T>&>;
    requires Numeric<T>;
    requires b.degree() >= 0;
};

// -----------------------------------------------------------------------------
// Curve concept
// -----------------------------------------------------------------------------

/**
 * CurveConcept
 *
 * A parametric curve in ℝ^n (typically ℝ^2 or ℝ^3).
 *
 * Requirements:
 *   - `evaluate(u)`  → NURBSPoint<T>  (control point evaluation at parameter u)
 *   - `degree()`    → int
 *   - `knot_vector()` → KnotVector<T>
 *   - `control_points()` → something range-like
 *   - `weights()`   → WeightVector<T>  (may be trivial uniform weights)
 *   - `parameter_domain()` → std::pair<T, T>
 *
 * The evaluate() function performs the De Boor or De Casteljau evaluation
 * and must return a point on the curve for any u in the parameter domain.
 */
template <typename C, typename T = double>
concept CurveConcept = ParametricEntity<C, T> && BasisFunctionProvider<C, T>
    && requires(const C& c, T u) {
        { c.evaluate(u) } -> std::same_as<NURBSPoint<T>>;
        { c.control_points() } -> std::same_as<std::vector<NURBSPoint<T>>>;
        { c.weights() } -> std::same_as<WeightVector<T>>;
};

// -----------------------------------------------------------------------------
// Surface concept
// -----------------------------------------------------------------------------

/**
 * SurfaceConcept
 *
 * A bivariate parametric surface in ℝ^n.
 *
 * Requirements:
 *   - `evaluate(u, v)` → NURBSPoint<T>
 *   - `u_degree()`     → int
 *   - `v_degree()`     → int
 *   - `u_knot_vector()` → KnotVector<T>
 *   - `v_knot_vector()` → KnotVector<T>
 *   - `control_points()` → 2D grid of NURBSPoint<T>
 *   - `weights()`      → WeightVector<T> (1D or 2D depending on representation)
 *   - `parameter_domain_u()` → std::pair<T, T>
 *   - `parameter_domain_v()` → std::pair<T, T>
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
    requires Numeric<T>;
    requires s.u_degree() >= 0;
    requires s.v_degree() >= 0;
};

// -----------------------------------------------------------------------------
// Point/vector arithmetic helpers
// -----------------------------------------------------------------------------

/// Concept for a geometric point or vector in ℝ^n
template <typename P, typename T = double>
concept PointLike = requires(const P& a, const P& b, T s) {
    { a + b } -> std::same_as<P>;
    { a - b } -> std::same_as<P>;
    { a * s } -> std::same_as<P>;
    { s * a } -> std::same_as<P>;
    { a / s } -> std::same_as<P>;
    requires Numeric<T>;
};

/// Homogeneous weight scalar companion to NURBSPoint
template <typename W, typename T = double>
concept WeightVector = requires(const W& w, T weight) {
    { w.size() } -> std::same_as<std::size_t>;
    { w[std::declval<std::size_t>()] } -> std::same_as<T>;
    requires Numeric<T>;
};

} // namespace nurbs::core

// numeric.hpp — Precision configuration and configurable tolerances
#pragma once

#include <limits>
#include <type_traits>
#include <cmath>

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Default tolerances — IEEE-754 double is assumed as the default scalar
// -----------------------------------------------------------------------------

/// Machine epsilon for float  (ISO IEEE-754 single)
inline constexpr float  default_epsilon_f32  = 1.1920929e-7f;

/// Machine epsilon for double (ISO IEEE-754 double)
inline constexpr double default_epsilon_f64  = 2.2204460492503131e-16;

/// Machine epsilon for long double (80-bit extended, typically)
inline constexpr long double default_epsilon_f80 = 1.0842021724855044e-19L;

/// Relative spacing between representable doubles (ulp(1))
inline constexpr double default_ulp_one = 2.0 * default_epsilon_f64;

// -----------------------------------------------------------------------------
// PrecisionConfig — scalar-type tag for compile-time dispatch
// -----------------------------------------------------------------------------

/**
 * PrecisionConfig — a tag class that bundles machine-epsilon and ulp(1) for
 * a given floating-point type.
 *
 * Use `PrecisionConfig<double>` (the default) when writing generic code that
 * must respect the scalar's native precision.
 *
 * Access via: `PrecisionConfig<T>::epsilon`, `PrecisionConfig<T>::ulp_one`.
 */
template <NumericScalar T>
struct PrecisionConfig {
    static constexpr T epsilon  = []() {
        if constexpr (std::same_as<T, float>)        return static_cast<T>(default_epsilon_f32);
        else if constexpr (std::same_as<T, double>)   return static_cast<T>(default_epsilon_f64);
        else                                           return static_cast<T>(default_epsilon_f80);
    }();

    static constexpr T ulp_one  = []() {
        if constexpr (std::same_as<T, float>)        return static_cast<T>(2.0f * default_epsilon_f32);
        else if constexpr (std::same_as<T, double>)   return static_cast<T>(2.0 * default_epsilon_f64);
        else                                           return static_cast<T>(2.0L * default_epsilon_f80);
    }();

    /// The "safe" tolerance for comparisons involving this scalar type.
    /// Equal to 10 × ulp(1) — catches two closest representable doubles
    /// while tolerating the normal rounding noise in B-spline evaluation.
    static constexpr T default_tolerance() noexcept { return static_cast<T>(10) * ulp_one; }
};

// -----------------------------------------------------------------------------
// Tolerance — value-level, configurable comparison tolerance
// -----------------------------------------------------------------------------

/**
 * Tolerance — a value-level tolerance wrapper that implements approximate
 * equality with a configurable relative and absolute component.
 *
 * Construction:
 *   Tolerance<double> tol;                  // default: 10 * ulp(1) ≈ 2.2e-15
 *   Tolerance<double> tol{1e-9};            // explicit absolute/relative cap
 *   Tolerance<double> tol{1e-6, 1e-12};     // rel cap, abs floor
 *
 * Usage:
 *   if (tol.eq(a, b))   // |a - b| ≤ rel * max(|a|,|b|) + abs
 *   if (tol.ne(a, b))   // negation of above
 *   if (tol.zero(x))    // |x| ≤ abs
 *   if (tol.gt(a, b))   // a > b + tolerance
 *   if (tol.lt(a, b))   // a < b - tolerance
 *   if (tol.ge(a, b))   // a ≥ b - tolerance
 *   if (tol.le(a, b))   // a ≤ b + tolerance
 */
template <NumericScalar T>
class Tolerance {
public:
    using scalar_type = T;

    /// Default tolerance: 10 × ulp(1) of the scalar type.
    static Tolerance defaults() noexcept {
        return Tolerance{PrecisionConfig<T>::default_tolerance()};
    }

    /// Construct with a single cap — used as both the relative cap and the
    /// absolute floor (the conservative choice when nothing is known about a/b).
    explicit Tolerance(T cap) : rel_cap_(cap), abs_floor_(cap) {}

    /// Construct with separate relative cap and absolute floor.
    Tolerance(T rel_cap, T abs_floor) : rel_cap_(rel_cap), abs_floor_(abs_floor) {}

    /// Construct from a PrecisionConfig tag.
    explicit Tolerance(PrecisionConfig<T>) noexcept
        : rel_cap_(PrecisionConfig<T>::default_tolerance())
        , abs_floor_(PrecisionConfig<T>::default_tolerance())
    {}

    // ----- comparison API -----

    /// Approximate equality: |a - b| ≤ rel_cap * max(|a|,|b|) + abs_floor_
    [[nodiscard]] bool eq(T a, T b) const noexcept {
        const T diff = a > b ? a - b : b - a;
        return diff <= rel_cap_ * (std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b))
                     + abs_floor_;
    }

    /// Approximate inequality (negation of eq).
    [[nodiscard]] bool ne(T a, T b) const noexcept { return !eq(a, b); }

    /// True when |x| ≤ abs_floor_  (near-zero test).
    [[nodiscard]] bool zero(T x) const noexcept {
        const T ax = x > T{0} ? x : -x;
        return ax <= abs_floor_;
    }

    /// a > b + tolerance
    [[nodiscard]] bool gt(T a, T b) const noexcept {
        const T diff = a - b;
        return diff > rel_cap_ * (std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b))
                    + abs_floor_;
    }

    /// a < b - tolerance
    [[nodiscard]] bool lt(T a, T b) const noexcept { return gt(b, a); }

    /// a ≥ b - tolerance
    [[nodiscard]] bool ge(T a, T b) const noexcept { return !lt(a, b); }

    /// a ≤ b + tolerance
    [[nodiscard]] bool le(T a, T b) const noexcept { return !gt(a, b); }

    // ----- accessors -----
    [[nodiscard]] T relative_cap()   const noexcept { return rel_cap_; }
    [[nodiscard]] T absolute_floor() const noexcept { return abs_floor_; }

private:
    T rel_cap_;   // relative tolerance multiplier
    T abs_floor_; // absolute tolerance floor
};

// -----------------------------------------------------------------------------
// Specialised standard-library overloads for Tolerance
// -----------------------------------------------------------------------------

/// Convenience alias matching the NURBS Book notation: `Toleranced` = tolerance for double.
using Toleranced = Tolerance<double>;

// Default Tolerance instance (double, 10 × ulp(1))
inline constexpr Toleranced default_tolerance = Toleranced::defaults();

// -----------------------------------------------------------------------------
// Numeric helpers
// -----------------------------------------------------------------------------

/// Safe division: returns 0 when divisor is near-zero under the given tolerance.
template <NumericScalar T>
[[nodiscard]] T safe_div(T a, T b, T zero_tol) noexcept {
    const T abs_b = b > T{0} ? b : -b;
    if (abs_b <= zero_tol) return T{0};
    return a / b;
}

/// Clamp x to [lo, hi]
template <NumericScalar T>
[[nodiscard]] T clamp(T x, T lo, T hi) noexcept {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/// Linear interpolation: a + t * (b - a),  t ∈ [0,1]
template <NumericScalar T>
[[nodiscard]] T lerp(T a, T b, T t) noexcept {
    return a + t * (b - a);
}

/// Maps u from [t0, t1] onto [s0, s1] (linear reparameterisation).
template <NumericScalar T>
[[nodiscard]] T remap(T u, T t0, T t1, T s0, T s1) noexcept {
    const T denom = t1 - t0;
    if (denom == T{0}) return s0; // degenerate domain → constant map
    return s0 + (u - t0) / denom * (s1 - s0);
}

// -----------------------------------------------------------------------------
// Compile-time safe comparisons
// -----------------------------------------------------------------------------

/// True when two floating-point literals are exactly equal at compile time.
template <NumericScalar T, T a, T b>
inline constexpr bool ct_eq_v = (a == b);

} // namespace nurbs::core

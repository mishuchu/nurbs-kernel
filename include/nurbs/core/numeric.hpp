// numeric.hpp — Precision configuration and configurable tolerances
#pragma once

#include <algorithm>
#include <cmath>
#include <limits>
#include <type_traits>
#include <utility>

namespace nurbs::core {

// -----------------------------------------------------------------------------
// NumericScalar_ — must be defined before PrecisionConfig and Tolerance use it
// -----------------------------------------------------------------------------

/**
 * NumericScalar_ — builtin floating-point types usable as NURBS coordinate
 * and parameter scalars.  Named NumericScalar_ here (not NumericScalar) so that
 * downstream headers can include this file and optionally re-export the name
 * without conflicts.
 */
template <typename T>
concept NumericScalar_ = std::floating_point<T>;

// -----------------------------------------------------------------------------
// Machine epsilon constants
// -----------------------------------------------------------------------------

inline constexpr double default_epsilon_f64  = 2.2204460492503131e-16;
inline constexpr float   default_epsilon_f32  = 1.1920929e-7f;
inline constexpr long double default_epsilon_f80 = 1.0842021724855044e-19L;

// -----------------------------------------------------------------------------
// PrecisionConfig — scalar-type tag for compile-time dispatch
// -----------------------------------------------------------------------------

/**
 * PrecisionConfig — a tag class that bundles machine-epsilon and ulp(1) for
 * a given floating-point type.
 *
 * Access via: `PrecisionConfig<T>::epsilon`, `PrecisionConfig<T>::ulp_one`.
 */
template <NumericScalar_ T>
struct PrecisionConfig {
    static constexpr T epsilon = []() {
        if constexpr (std::same_as<T, float>)             return static_cast<T>(default_epsilon_f32);
        else if constexpr (std::same_as<T, double>)       return static_cast<T>(default_epsilon_f64);
        else                                               return static_cast<T>(default_epsilon_f80);
    }();

    static constexpr T ulp_one = []() {
        if constexpr (std::same_as<T, float>)             return static_cast<T>(2.0f * default_epsilon_f32);
        else if constexpr (std::same_as<T, double>)       return static_cast<T>(2.0  * default_epsilon_f64);
        else                                               return static_cast<T>(2.0L  * default_epsilon_f80);
    }();

    /// "Safe" tolerance: 10 × ulp(1) — catches two closest representable
    /// floats while tolerating normal rounding noise in B-spline evaluation.
    static constexpr T default_tolerance() noexcept { return static_cast<T>(10) * ulp_one; }
};

// -----------------------------------------------------------------------------
// Tolerance — value-level, configurable comparison tolerance
// -----------------------------------------------------------------------------

/**
 * Tolerance — a value-level tolerance wrapper implementing approximate equality.
 *
 * Construction:
 *   Tolerance<double> tol;                     // default: 10 * ulp(1) ≈ 2.2e-15
 *   Tolerance<double> tol{1e-9};               // explicit cap
 *   Tolerance<double> tol{1e-6, 1e-12};       // rel_cap, abs_floor
 *
 * Usage:
 *   if (tol.eq(a, b))   // |a-b| ≤ rel_cap * max(|a|,|b|) + abs_floor
 *   if (tol.ne(a, b))   // negation
 *   if (tol.zero(x))    // |x| ≤ abs_floor
 *   if (tol.gt(a, b))   // a > b + tolerance
 *   if (tol.lt(a, b))   // a < b - tolerance
 *   if (tol.ge(a, b))   // a ≥ b - tolerance
 *   if (tol.le(a, b))   // a ≤ b + tolerance
 */
template <NumericScalar_ T>
class Tolerance {
public:
    using scalar_type = T;

    static constexpr Tolerance defaults() noexcept {
        return Tolerance{PrecisionConfig<T>::default_tolerance()};
    }

    explicit constexpr Tolerance(T cap) : rel_cap_(cap), abs_floor_(cap) {}

    constexpr Tolerance(T rel_cap, T abs_floor) : rel_cap_(rel_cap), abs_floor_(abs_floor) {}

    explicit constexpr Tolerance(PrecisionConfig<T>) noexcept
        : rel_cap_(PrecisionConfig<T>::default_tolerance())
        , abs_floor_(PrecisionConfig<T>::default_tolerance())
    {}

    // ----- comparison API -----

    /// |a - b| ≤ rel_cap * max(|a|,|b|) + abs_floor_
    [[nodiscard]] bool eq(T a, T b) const noexcept {
        const T diff = a > b ? a - b : b - a;
        return diff <= rel_cap_ * (std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b))
                     + abs_floor_;
    }

    [[nodiscard]] bool ne(T a, T b) const noexcept { return !eq(a, b); }

    [[nodiscard]] bool zero(T x) const noexcept {
        const T ax = x > T{0} ? x : -x;
        return ax <= abs_floor_;
    }

    [[nodiscard]] bool gt(T a, T b) const noexcept {
        const T diff = a - b;
        return diff > rel_cap_ * (std::abs(a) > std::abs(b) ? std::abs(a) : std::abs(b))
                    + abs_floor_;
    }

    [[nodiscard]] bool lt(T a, T b) const noexcept { return gt(b, a); }
    [[nodiscard]] bool ge(T a, T b) const noexcept { return !lt(a, b); }
    [[nodiscard]] bool le(T a, T b) const noexcept { return !gt(a, b); }

    // ----- accessors -----
    [[nodiscard]] T relative_cap()   const noexcept { return rel_cap_; }
    [[nodiscard]] T absolute_floor() const noexcept { return abs_floor_; }

private:
    T rel_cap_;   // relative tolerance multiplier
    T abs_floor_; // absolute tolerance floor
};

// -----------------------------------------------------------------------------
// Typedef aliases
// -----------------------------------------------------------------------------

using Toleranced = Tolerance<double>;

inline constexpr Toleranced default_tolerance = Toleranced::defaults();

// -----------------------------------------------------------------------------
// Numeric helpers
// -----------------------------------------------------------------------------

/// Safe division: returns 0 when divisor is near-zero under zero_tol.
template <NumericScalar_ T>
[[nodiscard]] T safe_div(T a, T b, T zero_tol) noexcept {
    const T abs_b = b > T{0} ? b : -b;
    if (abs_b <= zero_tol) return T{0};
    return a / b;
}

/// Clamp x to [lo, hi]
template <NumericScalar_ T>
[[nodiscard]] T clamp(T x, T lo, T hi) noexcept {
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

/// Linear interpolation: a + t * (b - a),  t ∈ [0,1]
template <NumericScalar_ T>
[[nodiscard]] T lerp(T a, T b, T t) noexcept {
    return a + t * (b - a);
}

/// Maps u from [t0, t1] onto [s0, s1] (linear reparameterisation).
template <NumericScalar_ T>
[[nodiscard]] T remap(T u, T t0, T t1, T s0, T s1) noexcept {
    const T denom = t1 - t0;
    if (denom == T{0}) return s0;
    return s0 + (u - t0) / denom * (s1 - s0);
}

} // namespace nurbs::core

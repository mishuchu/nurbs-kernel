// types.hpp — Core type aliases and classes for NURBS kernel
#pragma once

#include <array>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <initializer_list>

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Point — Cartesian point in ℝ^n
// -----------------------------------------------------------------------------

/**
 * Point — a Cartesian point in n-dimensional space.
 *
 * @tparam T      Floating-point scalar type (float, double, long double)
 * @tparam Dim    Compile-time dimension, or 0 for dynamic (runtime) dimension
 */
template <NumericScalar T = double, std::size_t Dim = 0>
class Point {
public:
    using scalar_type = T;

    static constexpr std::size_t static_dim = Dim;

    // ----- constructors -----
    constexpr Point() requires(Dim != 0) : coords_{} {}

    constexpr Point(std::initializer_list<T> coords_) : coords_{} {
        const auto n = coords_.size();
        if constexpr (Dim != 0) {
            if (n != Dim) throw std::invalid_argument(
                "Point: initializer list size " + std::to_string(n) +
                " does not match static dimension " + std::to_string(Dim));
        }
        std::copy(coords_.begin(), coords_.end(), coords_.begin());
    }

    // ----- dimension -----
    [[nodiscard]] constexpr std::size_t dimension() const noexcept {
        if constexpr (Dim != 0) return Dim;
        return coords_.size();
    }

    // ----- element access -----
    [[nodiscard]] T operator[](std::size_t i) const { return coords_[i]; }
    [[nodiscard]] T& operator[](std::size_t i) { return coords_[i]; }

    // ----- arithmetic -----
    [[nodiscard]] Point operator+(const Point& rhs) const {
        Point r = *this;
        for (std::size_t i = 0; i < r.dimension(); ++i) r.coords_[i] += rhs.coords_[i];
        return r;
    }

    [[nodiscard]] Point operator-(const Point& rhs) const {
        Point r = *this;
        for (std::size_t i = 0; i < r.dimension(); ++i) r.coords_[i] -= rhs.coords_[i];
        return r;
    }

    [[nodiscard]] Point operator*(T s) const {
        Point r = *this;
        for (std::size_t i = 0; i < r.dimension(); ++i) r.coords_[i] *= s;
        return r;
    }

    [[nodiscard]] Point operator/(T s) const {
        Point r = *this;
        for (std::size_t i = 0; i < r.dimension(); ++i) r.coords_[i] /= s;
        return r;
    }

    [[nodiscard]] Point operator-() const {
        Point r = *this;
        for (std::size_t i = 0; i < r.dimension(); ++i) r.coords_[i] = -r.coords_[i];
        return r;
    }

    // ----- euclidean norm -----
    [[nodiscard]] T norm() const {
        T sum = 0;
        for (std::size_t i = 0; i < dimension(); ++i) sum += coords_[i] * coords_[i];
        return std::sqrt(sum);
    }

    [[nodiscard]] T norm_squared() const {
        T sum = 0;
        for (std::size_t i = 0; i < dimension(); ++i) sum += coords_[i] * coords_[i];
        return sum;
    }

private:
    // Dynamic storage for Dim == 0, inline for Dim != 0
    std::vector<T> coords_{};
};

// template variable to distinguish Point from NURBSPoint
template <NumericScalar T, std::size_t N>
inline constexpr bool is_Point_v = false;

// -----------------------------------------------------------------------------
// NURBSPoint — Homogeneous-coordinate control point
// -----------------------------------------------------------------------------

/**
 * NURBSPoint — a control point in homogeneous form (Pw = w * P).
 *
 * In homogeneous coordinates a NURBS point is represented as (x, y, z, w)
 * where w is the weight. The Cartesian point is (x/w, y/w, z/w).
 *
 * This class stores both the weighted homogeneous coordinates and the
 * associated weight for efficiency during algorithm implementation.
 */
template <NumericScalar T = double>
class NURBSPoint {
public:
    using scalar_type = T;

    NURBSPoint() = default;

    /// Construct from weighted homogeneous coordinates and explicit weight.
    NURBSPoint(T x, T y, T w) : x_(x), y_(y), w_(w) {}

    /// Construct from weighted homogeneous coordinates with explicit weight for 3D.
    NURBSPoint(T x, T y, T z, T w) : x_(x), y_(y), z_(z), w_(w) {}

    /// Construct from a Cartesian Point and weight (auto-convert to homogeneous).
    explicit NURBSPoint(const Point<T>& p, T weight = T{1})
        : x_(p[0] * weight)
        , y_(p[1] * weight)
        , z_(p.size() > 2 ? p[2] * weight : T{0})
        , w_(weight)
    {}

    // ----- accessors -----
    [[nodiscard]] T x()  const noexcept { return x_; }
    [[nodiscard]] T y()  const noexcept { return y_; }
    [[nodiscard]] T z()  const noexcept { return z_; }
    [[nodiscard]] T w()  const noexcept { return w_; }  // weight

    /// Cartesian x = x/w
    [[nodiscard]] T cart_x() const noexcept { return w_ != T{0} ? x_ / w_ : x_; }
    /// Cartesian y = y/w
    [[nodiscard]] T cart_y() const noexcept { return w_ != T{0} ? y_ / w_ : y_; }
    /// Cartesian z = z/w
    [[nodiscard]] T cart_z() const noexcept { return w_ != T{0} ? z_ / w_ : z_; }

    // ----- arithmetic -----
    [[nodiscard]] NURBSPoint operator+(const NURBSPoint& rhs) const {
        return NURBSPoint(x_ + rhs.x_, y_ + rhs.y_, z_ + rhs.z_, w_ + rhs.w_);
    }

    [[nodiscard]] NURBSPoint operator-(const NURBSPoint& rhs) const {
        return NURBSPoint(x_ - rhs.x_, y_ - rhs.y_, z_ - rhs.z_, w_ - rhs.w_);
    }

    [[nodiscard]] NURBSPoint operator*(T s) const {
        return NURBSPoint(x_ * s, y_ * s, z_ * s, w_ * s);
    }

    [[nodiscard]] NURBSPoint operator/(T s) const {
        return NURBSPoint(x_ / s, y_ / s, z_ / s, w_ / s);
    }

    /// Dehomogenize → Cartesian Point (assumes 2D or 3D)
    [[nodiscard]] Point<T> cartesian(std::size_t dim = 2) const {
        if (dim == 2) return Point<T>{cart_x(), cart_y()};
        return Point<T>{cart_x(), cart_y(), cart_z()};
    }

private:
    T x_{0}, y_{0}, z_{0}, w_{1};
};

// -----------------------------------------------------------------------------
// KnotVector — Non-decreasing sequence of knot values
// -----------------------------------------------------------------------------

/**
 * KnotVector — a non-decreasing sequence of knot values.
 *
 * Represents the knot vector U = { u_0, u_1, …, u_m } used in B-spline / NURBS
 * curve and surface definitions.  The vector has m+1 knots for a basis of degree
 * p with n+1 control points where m = n + p + 1 (Algorithm A2.1 in NURBS Book).
 *
 * Invariant: knots are non-decreasing.
 */
template <NumericScalar T = double>
class KnotVector {
public:
    using value_type = T;
    using size_type  = std::size_t;

    // ----- constructors -----
    KnotVector() = default;

    explicit KnotVector(std::vector<T> knots) : knots_(std::move(knots)) {
        validate();
    }

    KnotVector(std::initializer_list<T> knots) : knots_(knots) {
        validate();
    }

    /// Uniform knot vector on [0,1] with `n+1` knots (clamped).
    static KnotVector uniform(std::size_t n_plus_1, T start = T{0}, T end = T{1}) {
        if (n_plus_1 < 2) throw std::invalid_argument("KnotVector::uniform: need at least 2 knots");
        std::vector<T> kv(n_plus_1);
        const T span = (end - start) / static_cast<T>(n_plus_1 - 1);
        for (std::size_t i = 0; i < n_plus_1; ++i) kv[i] = start + static_cast<T>(i) * span;
        return KnotVector(std::move(kv));
    }

    // ----- size -----
    [[nodiscard]] bool empty()    const noexcept { return knots_.empty(); }
    [[nodiscard]] size_type size() const noexcept { return knots_.size(); }

    // ----- access -----
    [[nodiscard]] const T& operator[](size_type i) const { return knots_[i]; }
    [[nodiscard]] const T& at(size_type i)        const { return knots_.at(i); }

    // ----- range helpers -----
    [[nodiscard]] const T& front() const noexcept { return knots_.front(); }
    [[nodiscard]] const T& back()  const noexcept { return knots_.back(); }

    [[nodiscard]] T  domain_start() const noexcept { return knots_.front(); }
    [[nodiscard]] T  domain_end()   const noexcept { return knots_.back(); }

    // ----- mutators -----
    void push_back(T knot) {
        if (!knots_.empty() && knot < knots_.back())
            throw std::invalid_argument("KnotVector: knot values must be non-decreasing");
        knots_.push_back(knot);
    }

    // ----- algorithm helpers -----
    /// Number of interior knots (knots strictly between domain boundaries)
    [[nodiscard]] size_type num_interior() const {
        return size() > 2 ? size() - 2 : 0;
    }

    /// Count multiplicity of knot at index i
    [[nodiscard]] size_type multiplicity(size_type i) const {
        size_type mult = 1;
        while (i + mult < size() && knots_[i + mult] == knots_[i]) ++mult;
        return mult;
    }

    /// Total multiplicity across the entire vector
    [[nodiscard]] size_type total_multiplicity() const {
        size_type total = 0;
        for (size_type i = 0; i < size(); ++i)
            total += multiplicity(i);
        return total;
    }

    // ----- iterators -----
    [[nodiscard]] auto begin()  const noexcept { return knots_.begin(); }
    [[nodiscard]] auto end()    const noexcept { return knots_.end(); }
    [[nodiscard]] auto begin()        noexcept { return knots_.begin(); }
    [[nodiscard]] auto end()          noexcept { return knots_.end(); }

private:
    void validate() const {
        for (std::size_t i = 1; i < knots_.size(); ++i) {
            if (knots_[i] < knots_[i - 1])
                throw std::invalid_argument(
                    "KnotVector: knot at index " + std::to_string(i) +
                    " is less than previous knot");
        }
    }

    std::vector<T> knots_;
};

// -----------------------------------------------------------------------------
// WeightVector — companion weights for NURBS control points
// -----------------------------------------------------------------------------

/**
 * WeightVector — a flat array of positive weights w_i associated with each
 * control point in a NURBS curve or surface.
 *
 * For rational curves the weights multiply the Cartesian control points to form
 * the homogeneous representation: Pw_i = w_i * P_i.
 * The WeightVector stores just the scalar weights; the product with the
 * control points is performed by NURBSPoint during construction.
 */
template <NumericScalar T = double>
class WeightVector {
public:
    using value_type = T;
    using size_type  = std::size_t;

    WeightVector() = default;

    explicit WeightVector(std::vector<T> weights) : weights_(std::move(weights)) {
        validate();
    }

    WeightVector(std::initializer_list<T> weights) : weights_(weights) {
        validate();
    }

    /// Uniform weights (all 1.0) for n control points
    static WeightVector uniform(size_type n) {
        return WeightVector(std::vector<T>(n, T{1}));
    }

    // ----- size -----
    [[nodiscard]] bool empty()    const noexcept { return weights_.empty(); }
    [[nodiscard]] size_type size() const noexcept { return weights_.size(); }

    // ----- access -----
    [[nodiscard]] const T& operator[](size_type i) const { return weights_[i]; }
    [[nodiscard]] const T& at(size_type i)        const { return weights_.at(i); }
    [[nodiscard]] T&       operator[](size_type i)       { return weights_[i]; }
    [[nodiscard]] T&       at(size_type i)              { return weights_.at(i); }

    // ----- mutators -----
    void push_back(T w) {
        if (w <= T{0})
            throw std::invalid_argument("WeightVector: weights must be positive");
        weights_.push_back(w);
    }

    // ----- algorithm helpers -----
    /// Sum of all weights
    [[nodiscard]] T sum() const {
        T s = T{0};
        for (auto w : weights_) s += w;
        return s;
    }

    /// Product of all weights
    [[nodiscard]] T product() const {
        T p = T{1};
        for (auto w : weights_) p *= w;
        return p;
    }

    // ----- iterators -----
    [[nodiscard]] auto begin() const noexcept { return weights_.begin(); }
    [[nodiscard]] auto end()   const noexcept { return weights_.end(); }

private:
    void validate() const {
        for (std::size_t i = 0; i < weights_.size(); ++i) {
            if (weights_[i] <= T{0})
                throw std::invalid_argument(
                    "WeightVector: weight at index " + std::to_string(i) +
                    " must be positive, got " + std::to_string(weights_[i]));
        }
    }

    std::vector<T> weights_;
};

// -----------------------------------------------------------------------------
// Type aliases for common use cases
// -----------------------------------------------------------------------------

/// 2D point (Cartesian)
template <NumericScalar T = double>
using Point2 = Point<T, 2>;

/// 3D point (Cartesian)
template <NumericScalar T = double>
using Point3 = Point<T, 3>;

/// 2D NURBS control point with weight
template <NumericScalar T = double>
using NURBSPoint2 = NURBSPoint<T>;

/// 3D NURBS control point with weight
template <NumericScalar T = double>
using NURBSPoint3 = NURBSPoint<T>;

// -----------------------------------------------------------------------------
// NumericScalar concept (used above)
// -----------------------------------------------------------------------------

/**
 * NumericScalar — builtin floating-point types that can be used as the
 * coordinate / parameter scalar for NURBS types.
 */
template <typename T>
concept NumericScalar = std::same_as<T, float> || std::same_as<T, double>
                      || std::same_as<T, long double>;

} // namespace nurbs::core

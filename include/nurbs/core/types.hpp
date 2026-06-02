// types.hpp — Core type aliases and classes for NURBS kernel
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "numeric.hpp"  // provides NumericScalar_ (concept), Tolerance, etc.

namespace nurbs::core {

// -----------------------------------------------------------------------------
// Point — Cartesian point in ℝ^n (static or dynamic dimension)
// -----------------------------------------------------------------------------

/**
 * Point — a Cartesian point in n-dimensional space.
 *
 * @tparam T   Floating-point scalar type (float, double, long double)
 * @tparam Dim Compile-time dimension; 0 = runtime (dynamic).
 */
template <NumericScalar_ T = double, std::size_t Dim = 0>
class Point {
public:
    using scalar_type = T;
    static constexpr std::size_t static_dim = Dim;

    // Default constructor: zero-initializes all coordinates.
    constexpr Point() : coords_{} {}

    // Initializer-list constructor: copies from the list into coords_.
    explicit Point(std::initializer_list<T> il) : coords_{} {
        const auto n = il.size();
        if constexpr (Dim != 0) {
            if (n != Dim)
                throw std::invalid_argument(
                    "Point: initializer size " + std::to_string(n)
                    + " != static dim " + std::to_string(Dim));
        }
        // Copy from il into coords_ using a simple loop (avoids iterator issues
        // with self-modifying containers in initialiser lists).
        auto it = il.begin();
        for (std::size_t i = 0; i < il.size(); ++i) coords_[i] = it[i];
    }

    [[nodiscard]] constexpr std::size_t dimension() const noexcept {
        if constexpr (Dim != 0) return Dim;
        return coords_.size();
    }

    [[nodiscard]] T operator[](std::size_t i) const { return coords_[i]; }
    [[nodiscard]] T& operator[](std::size_t i)       { return coords_[i]; }

    [[nodiscard]] Point operator+(const Point& rhs) const {
        Point r = *this;
        for (std::size_t i = 0; i < dimension(); ++i) r.coords_[i] += rhs.coords_[i];
        return r;
    }

    [[nodiscard]] Point operator-(const Point& rhs) const {
        Point r = *this;
        for (std::size_t i = 0; i < dimension(); ++i) r.coords_[i] -= rhs.coords_[i];
        return r;
    }

    [[nodiscard]] Point operator*(T s) const {
        Point r = *this;
        for (std::size_t i = 0; i < dimension(); ++i) r.coords_[i] *= s;
        return r;
    }

    [[nodiscard]] Point operator/(T s) const {
        Point r = *this;
        for (std::size_t i = 0; i < dimension(); ++i) r.coords_[i] /= s;
        return r;
    }

    [[nodiscard]] Point operator-() const {
        Point r = *this;
        for (std::size_t i = 0; i < dimension(); ++i) r.coords_[i] = -r.coords_[i];
        return r;
    }

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
    std::vector<T> coords_{};
};

// -----------------------------------------------------------------------------
// NURBSPoint — Homogeneous-coordinate control point Pw = w · P
// -----------------------------------------------------------------------------

/**
 * NURBSPoint — a control point in homogeneous form.
 *
 * Stored as (x, y, z, w) where w is the weight and (x/w, y/w, z/w) is the
 * Cartesian position.  The homogeneous representation simplifies the
 * De Boor / De Casteljau algorithms because affine combinations are linear.
 */
template <NumericScalar_ T = double>
class NURBSPoint {
public:
    using scalar_type = T;

    NURBSPoint() = default;

    NURBSPoint(T x, T y, T w) : x_(x), y_(y), z_{}, w_(w) {}

    NURBSPoint(T x, T y, T z, T w) : x_(x), y_(y), z_(z), w_(w) {}

    /** Build from a Cartesian point and weight (converts to homogeneous). */
    template <std::size_t D>
    explicit NURBSPoint(const Point<T, D>& p, T weight = T{1}) {
        x_ = p[0] * weight;
        y_ = (D >= 2 && p.dimension() > 1) ? p[1] * weight : T{0};
        z_ = (D >= 3 && p.dimension() > 2) ? p[2] * weight : T{0};
        w_ = weight;
    }

    // ----- accessors -----
    [[nodiscard]] T x()  const noexcept { return x_; }
    [[nodiscard]] T y()  const noexcept { return y_; }
    [[nodiscard]] T z()  const noexcept { return z_; }
    [[nodiscard]] T w()  const noexcept { return w_; }

    [[nodiscard]] T cart_x() const noexcept { return w_ != T{0} ? x_ / w_ : x_; }
    [[nodiscard]] T cart_y() const noexcept { return w_ != T{0} ? y_ / w_ : y_; }
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

    /** Dehomogenize → Cartesian Point. */
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
 * KnotVector — a non-decreasing sequence of knot values U = {u_0 … u_m}.
 *
 * For a B-spline of degree p with n+1 control points the knot vector has
 * m+1 = n+p+2 knots (Algorithm A2.1, NURBS Book).  Clamped curves have
 * p+1 equal knots at each end of the domain.
 */
template <NumericScalar_ T = double>
class KnotVector {
public:
    using value_type = T;
    using size_type  = std::size_t;

    KnotVector() = default;

    explicit KnotVector(std::vector<T> knots) : knots_(std::move(knots)) { validate(); }
    explicit KnotVector(std::initializer_list<T> knots) : knots_(knots) { validate(); }

    /** Uniform knot vector on [start, end] with n+1 equally-spaced knots. */
    static KnotVector uniform(size_type n_plus_1, T start = T{0}, T end = T{1}) {
        if (n_plus_1 < 2) throw std::invalid_argument("KnotVector::uniform: need ≥ 2 knots");
        std::vector<T> kv(n_plus_1);
        const T span = (end - start) / static_cast<T>(n_plus_1 - 1);
        for (size_type i = 0; i < n_plus_1; ++i)
            kv[i] = start + static_cast<T>(i) * span;
        return KnotVector(std::move(kv));
    }

    // ----- size -----
    [[nodiscard]] bool empty()    const noexcept { return knots_.empty(); }
    [[nodiscard]] size_type size() const noexcept { return knots_.size(); }

    // ----- access -----
    [[nodiscard]] const T& operator[](size_type i) const { return knots_[i]; }
    [[nodiscard]] const T& at(size_type i)        const { return knots_.at(i); }

    [[nodiscard]] const T& front() const noexcept { return knots_.front(); }
    [[nodiscard]] const T& back()  const noexcept { return knots_.back(); }

    [[nodiscard]] T domain_start() const noexcept { return knots_.front(); }
    [[nodiscard]] T domain_end()   const noexcept { return knots_.back(); }

    // ----- mutators -----
    void push_back(T knot) {
        if (!knots_.empty() && knot < knots_.back())
            throw std::invalid_argument("KnotVector: knots must be non-decreasing");
        knots_.push_back(knot);
    }

    // ----- helpers -----
    /** Number of interior knots (strictly between domain boundaries). */
    [[nodiscard]] size_type num_interior() const {
        return size() > 2 ? size() - 2 : 0;
    }

    /** Multiplicity of knot at index i. */
    [[nodiscard]] size_type multiplicity(size_type i) const {
        size_type m = 1;
        while (i + m < size() && knots_[i + m] == knots_[i]) ++m;
        return m;
    }

    /** Total multiplicity across all knots. */
    [[nodiscard]] size_type total_multiplicity() const {
        size_type total = 0;
        for (size_type i = 0; i < size(); ++i) total += multiplicity(i);
        return total;
    }

    // ----- iterators -----
    [[nodiscard]] auto begin() const noexcept { return knots_.begin(); }
    [[nodiscard]] auto end()   const noexcept { return knots_.end(); }

private:
    void validate() const {
        for (size_type i = 1; i < knots_.size(); ++i) {
            if (knots_[i] < knots_[i - 1])
                throw std::invalid_argument(
                    "KnotVector: knot[" + std::to_string(i) + "] < knot["
                    + std::to_string(i - 1) + "]");
        }
    }

    std::vector<T> knots_;
};

// -----------------------------------------------------------------------------
// WeightVector — Companion weights for NURBS control points
// -----------------------------------------------------------------------------

/**
 * WeightVector — a flat array of positive weights w_i associated with each
 * control point in a NURBS curve or surface.
 *
 * The weights multiply the Cartesian control points to form the homogeneous
 * representation: Pw_i = w_i · P_i.
 */
template <NumericScalar_ T = double>
class WeightVector {
public:
    using value_type = T;
    using size_type  = std::size_t;

    WeightVector() = default;

    explicit WeightVector(std::vector<T> weights) : weights_(std::move(weights)) { validate(); }
    explicit WeightVector(std::initializer_list<T> weights) : weights_(weights) { validate(); }

    /** Uniform weights (all 1.0) for n control points. */
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
    [[nodiscard]] T&       at(size_type i)               { return weights_.at(i); }

    // ----- mutators -----
    void push_back(T w) {
        if (w <= T{0})
            throw std::invalid_argument(
                "WeightVector: weight must be positive, got " + std::to_string(w));
        weights_.push_back(w);
    }

    // ----- helpers -----
    [[nodiscard]] T sum()     const { T s{}; for (auto w : weights_) s += w; return s; }
    [[nodiscard]] T product() const { T p{T{1}}; for (auto w : weights_) p *= w; return p; }

    // ----- iterators -----
    [[nodiscard]] auto begin() const noexcept { return weights_.begin(); }
    [[nodiscard]] auto end()   const noexcept { return weights_.end(); }

private:
    void validate() const {
        for (size_type i = 0; i < weights_.size(); ++i) {
            if (weights_[i] <= T{0})
                throw std::invalid_argument(
                    "WeightVector: weights must be positive, got "
                    + std::to_string(weights_[i]) + " at index " + std::to_string(i));
        }
    }

    std::vector<T> weights_;
};

// -----------------------------------------------------------------------------
// Common type aliases
// -----------------------------------------------------------------------------

template <NumericScalar_ T = double>
using Point2 = Point<T, 2>;

template <NumericScalar_ T = double>
using Point3 = Point<T, 3>;

template <NumericScalar_ T = double>
using NURBSPoint2 = NURBSPoint<T>;

template <NumericScalar_ T = double>
using NURBSPoint3 = NURBSPoint<T>;

} // namespace nurbs::core

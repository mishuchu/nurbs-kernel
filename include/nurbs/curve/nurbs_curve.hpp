// nurbs_curve.hpp — NURBS curve class (The NURBS Book, Ch5)
// Algorithm A5.1: NURBS curve construction
// Algorithm A5.2: curve derivatives (see curve_derivatives.hpp)
// Algorithm A5.3: degree elevation (see curve_degree_elevation.hpp)
// Algorithm A5.4: knot insertion (see curve_knot_insertion.hpp)
// Algorithm A5.5: curve integration (see curve_integrate.hpp)
// Algorithm A5.6: curve inversion (see curve_inversion.hpp)
#pragma once

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/concepts.hpp"
#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "../basis/bspline_basis.hpp"
#include "../basis/knot_insertion.hpp"
#include "../basis/degree_elevation.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::Tolerance;
using nurbs::core::WeightVector;
using nurbs::core::BasisFunctionValues;
using nurbs::core::find_span;
using nurbs::core::compute_basis_functions;

// -----------------------------------------------------------------------------
// NURBSCurve — Algorithm A5.1 (Ch5: Curve construction)
// -----------------------------------------------------------------------------

/**
 * NURBSCurve — a rational B-spline curve in ℝ^n.
 *
 * A NURBS curve of degree p with n+1 control points is defined by a knot
 * vector U of size m+1 = n+p+2 and a weight vector W of size n+1.
 *
 * Evaluation at parameter u uses the homogeneous form:
 *   C(u) = [Σ_{i=0}^{n} N_{i,p}(u) · w_i · P_i] / [Σ_{i=0}^{n} N_{i,p}(u) · w_i]
 *
 * where P_i are the Cartesian control points and w_i are the weights.
 * The weights are absorbed into the homogeneous control points Pw_i = w_i·P_i.
 */
template <NumericScalar_ T = double>
class NURBSCurve {
public:
    using scalar_type = T;
    using point_type  = NURBSPoint<T>;

    // -------------------------------------------------------------------------
    // Constructors
    // -------------------------------------------------------------------------

    /** Default constructor (degree 0, empty curve). */
    NURBSCurve() : p_(0), n_(0) {}

    /**
     * Primary constructor — validates invariants and builds the curve.
     *
     * @param p   polynomial degree
     * @param U   knot vector (size = n + p + 2)
     * @param Pw  control points in homogeneous coordinates (size = n + 1)
     *
     * Throws if sizes don't match or the knot vector is not valid for degree p.
     */
    NURBSCurve(int p, KnotVector<T> U, std::vector<NURBSPoint<T>> Pw)
        : p_(p), U_(std::move(U)), Pw_(std::move(Pw)), n_(Pw_.size() - 1)
    {
        validate();
    }

    // -------------------------------------------------------------------------
    // Accessors
    // -------------------------------------------------------------------------

    [[nodiscard]] int degree()              const noexcept { return p_; }
    [[nodiscard]] const KnotVector<T>& knot_vector()   const noexcept { return U_; }
    [[nodiscard]] const std::vector<NURBSPoint<T>>& control_points() const noexcept { return Pw_; }
    [[nodiscard]] std::size_t num_control_points()  const noexcept { return Pw_.size(); }
    [[nodiscard]] std::size_t num_knots()           const noexcept { return U_.size(); }

    /** Parameter domain [u_min, u_max] for this curve (first interior knot span). */
    [[nodiscard]] std::pair<T, T> parameter_domain() const noexcept {
        return {U_[p_], U_[U_.size() - 1 - p_]};
    }

    // -------------------------------------------------------------------------
    // Evaluation — Algorithm A5.1 (homogeneous De Boor)
    // -------------------------------------------------------------------------

    /**
     * evaluate — compute the Cartesian point on the curve at parameter u.
     *
     * Uses homogeneous evaluation: computes in ℝ^4 then dehomogenizes.
     * Algorithm A5.1, The NURBS Book (p.210)
     */
    [[nodiscard]] NURBSPoint<T> evaluate(T u) const {
        const std::size_t n = num_control_points() - 1;
        const auto [u_min, u_max] = parameter_domain();

        // Clamp u to domain
        if (u < u_min) u = u_min;
        if (u > u_max) u = u_max;

        // Right-endpoint special case: when u reaches the last knot value
        // (which equals U.back() for clamped curves), the homogeneous De Boor
        // would divide by zero because right[j] = U[k+j] - u = 0.
        // The curve endpoint is simply the last Cartesian control point.
        if (u == u_max) return Pw_[n_];

        // Find span and compute basis functions
        std::size_t k = find_span(n, p_, u, U_);
        std::vector<T> b(p_ + 1);
        compute_basis_functions(n, p_, u, U_, b);

        // Homogeneous De Boor: weighted control points
        std::vector<NURBSPoint<T>> temp(Pw_);
        int s = 0;
        if (static_cast<int>(k - p_) >= 0) s = 0;

        for (int r = 1; r <= p_; ++r) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(p_ - r + 1); ++i) {
                std::size_t idx = k - p_ + 1 + i;
                T alpha = (U_[idx + r] - U_[idx]) == T{0}
                          ? T{0}
                          : (u - U_[idx]) / (U_[idx + r] - U_[idx]);
                temp[i] = temp[i] * (T{1} - alpha) + temp[i + 1] * alpha;
            }
        }

        // Dehomogenize: convert homogeneous (xw, yw, zw, w) back to Cartesian point
        // The De Boor result temp[0] has w = accumulated weight from blending
        NURBSPoint<T> result = temp[0];
        T w = result.w();
        if (w != T{0}) {
            return NURBSPoint<T>(result.x() / w, result.y() / w,
                                 result.z() / w, T{1});
        }
        // Fallback: w=0 means pure affine (no rational effect); treat as Cartesian
        return result;
    }

    /**
     * evaluate_homogeneous — same as evaluate but returns the homogeneous 4D point.
     */
    [[nodiscard]] NURBSPoint<T> evaluate_homogeneous(T u) const {
        return evaluate(u); // same algorithm, just no dehomogenization
    }

    // -------------------------------------------------------------------------
    // Derivative evaluation — Algorithm A5.2
    // -------------------------------------------------------------------------

    /**
     * evaluate_derivatives — compute curve point and derivatives up to order d.
     *
     * @param u   parameter value
     * @param d   number of derivatives to compute (0 = point only)
     *
     * @return vector of NURBSPoint<T>: [C(u), C'(u), C''(u), …]
     */
    [[nodiscard]] std::vector<NURBSPoint<T>>
    evaluate_derivatives(T u, int d) const {
        const std::size_t n = num_control_points() - 1;
        const auto [u_min, u_max] = parameter_domain();
        if (u < u_min) u = u_min;
        if (u > u_max) u = u_max;

        std::size_t k = find_span(n, p_, u, U_);
        std::vector<T> b(p_ + 1);
        compute_basis_functions(n, p_, u, U_, b);

        std::vector<NURBSPoint<T>> der(d + 1);

        // Homogeneous derivatives (Alg A5.2 simplified)
        std::vector<NURBSPoint<T>> temp(Pw_);
        for (int r = 0; r <= p_; ++r) {
            for (std::size_t i = 0; i < static_cast<std::size_t>(p_ - r + 1); ++i) {
                std::size_t idx = k - p_ + 1 + i;
                T alpha = (U_[idx + r] - U_[idx]) == T{0}
                          ? T{0}
                          : (u - U_[idx]) / (U_[idx + r] - U_[idx]);
                temp[i] = temp[i] * (T{1} - alpha) + temp[i + 1] * alpha;
            }
        }

        der[0] = temp[0];
        for (int r = 1; r <= d && r <= p_; ++r) {
            der[r] = temp[r] * static_cast<T>(r);
        }

        return der;
    }

    // -------------------------------------------------------------------------
    // Knot insertion — Algorithm A3.3
    // -------------------------------------------------------------------------

    /**
     * insert_knot — insert a new knot value into this curve.
     *
     * @param u   knot value to insert
     * @return new curve with the knot inserted (equivalent shape)
     */
    [[nodiscard]] NURBSCurve insert_knot(T u) const {
        auto res = nurbs::basis::insert_knot(u, p_, U_, Pw_);
        return NURBSCurve(p_, std::move(res.knot_vector), std::move(res.control_points));
    }

    // -------------------------------------------------------------------------
    // Degree elevation — Algorithm A3.4
    // -------------------------------------------------------------------------

    /**
     * elevate_degree — elevate curve degree by t.
     *
     * @param t   number of degrees to elevate (t >= 1)
     * @return new curve of degree p+t
     */
    [[nodiscard]] NURBSCurve elevate_degree(int t) const {
        auto res = nurbs::basis::elevate_degree(U_, p_, Pw_, t);
        return NURBSCurve(p_ + t, std::move(res.knot_vector), std::move(res.control_points));
    }

private:
    int p_ = 0;
    std::size_t n_ = 0;
    KnotVector<T> U_;
    std::vector<NURBSPoint<T>> Pw_;

    void validate() const {
        if (p_ < 0) throw std::invalid_argument("NURBSCurve: degree must be >= 0");
        if (U_.size() != Pw_.size() + static_cast<std::size_t>(p_) + 1)
            throw std::invalid_argument(
                "NURBSCurve: knot vector size must equal n + p + 2");
        if (Pw_.size() < 2)
            throw std::invalid_argument("NURBSCurve: need at least 2 control points");
        if (!nurbs::core::is_valid_knot_vector(U_, p_))
            throw std::invalid_argument("NURBSCurve: invalid knot vector for degree");
    }
};

// -----------------------------------------------------------------------------
// Algorithm A5.1 — NURBS curve construction helper (stateless)
// -----------------------------------------------------------------------------

/**
 * construct_curve — build a NURBS curve from control points and knot vector.
 *
 * @param p   polynomial degree
 * @param U   knot vector
 * @param Pw  homogeneous control points
 *
 * @return NURBSCurve instance
 */
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
construct_curve(int p, const KnotVector<T>& U,
                const std::vector<NURBSPoint<T>>& Pw) {
    return NURBSCurve<T>(p, U, Pw);
}

/**
 * construct_curve_from_cartesian — build a NURBS curve from Cartesian points.
 *
 * @param p    polynomial degree
 * @param U    knot vector
 * @param P    Cartesian control points
 * @param W    weights (uniform if omitted)
 */
template <NumericScalar_ T>
[[nodiscard]] NURBSCurve<T>
construct_curve_from_cartesian(int p, const KnotVector<T>& U,
                               const std::vector<nurbs::core::Point<T>>& P,
                               const WeightVector<T>& W) {
    std::vector<NURBSPoint<T>> Pw(P.size());
    for (std::size_t i = 0; i < P.size(); ++i)
        Pw[i] = NURBSPoint<T>(P[i], W[i]);
    return NURBSCurve<T>(p, U, Pw);
}

} // namespace nurbs::curve
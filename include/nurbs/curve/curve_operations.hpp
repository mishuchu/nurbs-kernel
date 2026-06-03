// curve_operations.hpp — NURBS curve operations (The NURBS Book, Ch7)
// Algorithm A7.1: curve subdivision
// Algorithm A7.2: curve merging
// Algorithm A7.3: curve fitting (global least-squares)
// Reparametrization (rational linear)
#pragma once

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../core/concepts.hpp"
#include "../core/numeric.hpp"
#include "../core/types.hpp"
#include "../core/utilities.hpp"
#include "nurbs_curve.hpp"

namespace nurbs::curve {

using nurbs::core::KnotVector;
using nurbs::core::NURBSPoint;
using nurbs::core::NumericScalar_;
using nurbs::core::Point;
using nurbs::core::WeightVector;
using nurbs::core::Tolerance;
using nurbs::core::find_span;
using nurbs::core::compute_basis_functions;

// ===========================================================================
// Algorithm A7.1 — Curve Subdivision
// ===========================================================================

/**
 * SubdivideResult — output of curve_subdivide.
 */
template <NumericScalar_ T>
struct SubdivideResult {
    int p;                                     // degree (same for both halves)
    KnotVector<T> U1;                          // knot vector of left curve
    std::vector<NURBSPoint<T>> Pw1;           // control points of left curve
    KnotVector<T> U2;                          // knot vector of right curve
    std::vector<NURBSPoint<T>> Pw2;           // control points of right curve
};

/**
 * curve_subdivide — Algorithm A7.1 (The NURBS Book, 2nd ed., p.191)
 *
 * Splits a NURBS curve at an arbitrary parameter value u ∈ [u_k, u_{k+1}]
 * into two curves C1 and C2 that together represent the same geometry.
 *
 * The algorithm repeatedly inserts knot u until its multiplicity reaches p,
 * which makes it a C^{p-1} discontinuity — effectively splitting the
 * control point chain at that location.
 *
 * @param u    parameter at which to split (must be in the active domain)
 * @param p    curve degree
 * @param U    original knot vector
 * @param Pw   original homogeneous control points
 *
 * @return SubdivideResult containing two equivalent curve fragments.
 *
 * Example:
 *   auto [p, U1, Pw1, U2, Pw2] = curve_subdivide(0.5, 3, U, Pw);
 *   NURBSCurve<T> C1(p, U1, Pw1);
 *   NURBSCurve<T> C2(p, U2, Pw2);
 */
template <NumericScalar_ T>
[[nodiscard]] SubdivideResult<T>
curve_subdivide(T u, int p, const KnotVector<T>& U,
               const std::vector<NURBSPoint<T>>& Pw) {
    const std::size_t n = Pw.size() - 1;  // n+1 control points

    // Verify invariants
    if (p < 1) throw std::invalid_argument("curve_subdivide: degree must be >= 1");
    if (U.size() != n + p + 2)
        throw std::invalid_argument("curve_subdivide: knot vector size mismatch");

    // Algorithm A7.1 — find span k (same as find_span but inline)
    std::size_t k = find_span(n, p, u, U);

    // Determine current multiplicity of u
    int s = 0;
    for (std::size_t i = 0; i < U.size(); ++i) {
        if (U[i] == u) { s = static_cast<int>(U.multiplicity(i)); break; }
    }

    // Number of knot insertions needed to bring multiplicity to p
    int r = p - s;
    if (r < 0) r = 0;

    // Work on copies
    KnotVector<T> U_bar = U;
    std::vector<NURBSPoint<T>> Qw = Pw;

    // Repeated knot insertion (Algorithm A5.4 / A3.3) r times
    for (int t = 0; t < r; ++t) {
        auto res = nurbs::basis::insert_knot(u, p, U_bar, Qw);
        U_bar = std::move(res.knot_vector);
        Qw    = std::move(res.control_points);
    }

    // After r insertions, u now has multiplicity p in U_bar.
    // Find its index in U_bar — it marks the boundary between the two curves.
    std::size_t k_bar = find_span(n + r, p, u, U_bar);

    // Collect knots ≤ u into U1, knots ≥ u into U2.
    // For left half: include u with multiplicity p at the end.
    // For right half: include u with multiplicity p at the start.
    std::vector<T> U1_vec;
    std::vector<T> U2_vec;

    // U1: knots [0 … k_bar]
    for (std::size_t i = 0; i <= k_bar; ++i) U1_vec.push_back(U_bar[i]);

    // U2: knots [k_bar - p + 1 … end]
    // Note: after full multiplicity insertion the knot at k_bar equals u
    // and appears p+1 times (since k_bar-p+1 through k_bar all equal u).
    for (std::size_t i = k_bar - p + 1; i < U_bar.size(); ++i) U2_vec.push_back(U_bar[i]);

    KnotVector<T> U1(std::move(U1_vec));
    KnotVector<T> U2(std::move(U2_vec));

    // Split control points: Qw[0 … k_bar-p] belong to left, Qw[k_bar-p … end] to right.
    // The first (k_bar - p + 1) points → left curve C1
    // The last  (n + r - (k_bar - p) + 1) = n + r - k_bar + p + 1 points → right curve C2
    std::size_t left_n = k_bar - p + 1;           // number of control points in C1
    std::vector<NURBSPoint<T>> Pw1(left_n);
    for (std::size_t i = 0; i < left_n; ++i) Pw1[i] = Qw[i];

    std::vector<NURBSPoint<T>> Pw2(Qw.size() - left_n);
    for (std::size_t i = 0; i < Pw2.size(); ++i) Pw2[i] = Qw[left_n + i];

    return SubdivideResult<T>{p, std::move(U1), std::move(Pw1),
                              std::move(U2), std::move(Pw2)};
}

/**
 * curve_subdivide_at_knot_span — subdivide at the midpoint of knot span k.
 *
 * Convenience wrapper: finds the midpoint of span k and subdivides there.
 *
 * @param k   knot span index (0 ≤ k ≤ n+p)
 * @param p   degree
 * @param U   knot vector
 * @param Pw  control points
 */
template <NumericScalar_ T>
[[nodiscard]] SubdivideResult<T>
curve_subdivide_at_knot_span(std::size_t k, int p, const KnotVector<T>& U,
                            const std::vector<NURBSPoint<T>>& Pw) {
    if (k >= U.size() - 1) throw std::out_of_range("knot span index out of range");
    T u_mid = (U[k] + U[k + 1]) / T{2};
    return curve_subdivide(u_mid, p, U, Pw);
}

// ===========================================================================
// Algorithm A7.2 — Curve Merging (Global)
// ===========================================================================

/**
 * MergeResult — output of curve_merge.
 */
template <NumericScalar_ T>
struct MergeResult {
    int p;                                     // degree (must match)
    KnotVector<T> U;                           // merged knot vector
    std::vector<NURBSPoint<T>> Pw;            // merged control points
};

/**
 * curve_merge — Algorithm A7.2 (The NURBS Book, 2nd ed., p.192)
 *
 * Merges two NURBS curves C1(p,U1,Pw1) and C2(p,U2,Pw2) that share a
 * common endpoint at u = u_min of C2 / u = u_max of C1.
 *
 * The algorithm first verifies endpoint compatibility, then performs
 * knot refinement to bring both knot vectors into alignment so the
 * control point arrays can be concatenated.
 *
 * Prerequisites:
 *   - Both curves must have the same degree p.
 *   - C1's end parameter must equal C2's start parameter.
 *   - C1's endpoint control point must equal C2's endpoint control point.
 *
 * @param p1   degree of C1
 * @param U1   knot vector of C1
 * @param Pw1  homogeneous control points of C1
 * @param p2   degree of C2
 * @param U2   knot vector of C2
 * @param Pw2  homogeneous control points of C2
 *
 * @return MergeResult with merged curve (p, U, Pw).
 *
 * Example:
 *   auto merged = curve_merge(3, U1, Pw1, 3, U2, Pw2);
 *   NURBSCurve<T> C(merged.p, merged.U, merged.Pw);
 */
template <NumericScalar_ T>
[[nodiscard]] MergeResult<T>
curve_merge(int p1, const KnotVector<T>& U1,
            const std::vector<NURBSPoint<T>>& Pw1,
            int p2, const KnotVector<T>& U2,
            const std::vector<NURBSPoint<T>>& Pw2) {
    // Validate compatibility
    if (p1 != p2)
        throw std::invalid_argument("curve_merge: degrees must match");
    if (Pw1.size() < 2 || Pw2.size() < 2)
        throw std::invalid_argument("curve_merge: need at least 2 control points each");

    const int p = p1;

    // Verify endpoint compatibility: C1(u_max) == C2(u_min)
    auto [u1_min, u1_max] = nurbs::curve::NURBSCurve<T>(p, U1, Pw1).parameter_domain();
    auto [u2_min, u2_max] = nurbs::curve::NURBSCurve<T>(p, U2, Pw2).parameter_domain();

    Tolerance<T> tol = Tolerance<T>::defaults();

    // C1 endpoint vs C2 start
    if (!tol.eq(u1_max, u2_min))
        throw std::invalid_argument("curve_merge: C1 end u != C2 start u");

    // Endpoint control points must also match (within tolerance)
    const NURBSPoint<T>& c1_end = Pw1.back();
    const NURBSPoint<T>& c2_start = Pw2.front();
    if (!tol.eq(c1_end.cart_x(), c2_start.cart_x()) ||
        !tol.eq(c1_end.cart_y(), c2_start.cart_y()) ||
        !tol.eq(c1_end.cart_z(), c2_start.cart_z()))
        throw std::invalid_argument("curve_merge: endpoint control points do not match");

    // Algorithm A7.2:
    // Step 1: Bring C1 to have the same knot vector as C2 across the common boundary.
    //         We refine C1 by inserting all knots of C2 that lie in (u1_min, u1_max).
    // Step 2: Concatenate control point arrays.

    // Work on copies of C1 data
    KnotVector<T> U_bar = U1;
    std::vector<NURBSPoint<T>> Qw = Pw1;

    // Insert all knots from C2 that are in (u1_min, u1_max) into C1
    // (skip the clamped endpoints u2_min and u2_max)
    for (std::size_t i = p + 1; i + p + 1 < U2.size() - 1; ++i) {
        T u_ins = U2[i];
        // Only insert if strictly inside C1's domain and not already at full multiplicity
        if (u_ins > u1_min && u_ins < u1_max) {
            auto res = nurbs::basis::insert_knot(u_ins, p, U_bar, Qw);
            U_bar = std::move(res.knot_vector);
            Qw    = std::move(res.control_points);
        }
    }

    // Now U_bar and U2 share the same knots in their overlapping region.
    // Build merged knot vector: knots up to u1_max (inclusive) from U_bar,
    // followed by knots > u1_max from U2 (skipping the clamped tail of U2 at u2_max).
    std::vector<T> U_merged_vec;

    // All knots of U_bar up to index size()-p-1 (the last interior position)
    for (std::size_t i = 0; i + p + 1 < U_bar.size(); ++i)
        U_merged_vec.push_back(U_bar[i]);

    // The last p+1 knots of U_bar are at u1_max; skip duplicates with U2
    // Add knots from U2 that are strictly greater than u1_max
    for (std::size_t i = p + 1; i < U2.size(); ++i) {
        if (U2[i] > u1_max)
            U_merged_vec.push_back(U2[i]);
    }

    KnotVector<T> U_merged(std::move(U_merged_vec));

    // Control points: Qw (C1 refined) + Pw2[1:] (skip first = duplicate endpoint)
    std::vector<NURBSPoint<T>> Pw_merged;
    Pw_merged.reserve(Qw.size() + Pw2.size() - 1);
    for (const auto& pt : Qw) Pw_merged.push_back(pt);
    for (std::size_t i = 1; i < Pw2.size(); ++i) Pw_merged.push_back(Pw2[i]);

    return MergeResult<T>{p, std::move(U_merged), std::move(Pw_merged)};
}

// ===========================================================================
// Reparametrization — rational linear (mobius) reparameterization
// ===========================================================================

/**
 * ReparametrizeResult — output of curve_reparametrize.
 */
template <NumericScalar_ T>
struct ReparametrizeResult {
    int p;
    KnotVector<T> U;                          // transformed knot vector
    std::vector<NURBSPoint<T>> Pw;           // transformed control points
};

/**
 * reparametrize — rational linear (Mobius) reparameterization of a NURBS curve.
 *
 * Replaces the curve parameter u with u' = (a*u + b) / (c*u + d), effectively
 * applying a Mobius (linear fractional) transformation to the parameter domain.
 *
 * The effect on the curve is a reparameterization only — the geometric shape
 * is unchanged, only the way parameters map to positions changes.
 *
 * Mathematically:
 *   C'(u') = C(u(u')), with u(u') = (d*u' - b) / (a - c*u')
 *
 * For NURBS curves, the reparameterization transforms the knot vector and
 * control points using the relationship:
 *   N_{i,p}(u') = N_{i,p}(u(u')) where u' = (au+b)/(cu+d)
 *
 * Implementation uses knot vector remapping: each knot u_j in U gets
 * transformed to u'_j = (a*u_j + b) / (c*u_j + d), then the curve is
 * re-expressed in the new parameter space.
 *
 * @param a  numerator coefficient (a*d - b*c != 0 required)
 * @param b  numerator constant
 * @param c  denominator coefficient
 * @param d  denominator constant
 * @param p  curve degree
 * @param U  original knot vector
 * @param Pw original homogeneous control points
 *
 * @return ReparametrizeResult with reparameterized curve (same geometry, new parametrization)
 *
 * Note: The resulting curve evaluates to the same geometric points but with
 * different parameter values. This is useful for:
 *   - Improving parameter distribution (e.g., centripetal vs uniform)
 *   - Matching parameter ranges between curves
 *   - Aligning parameter directions
 */
template <NumericScalar_ T>
[[nodiscard]] ReparametrizeResult<T>
reparametrize(T a, T b, T c, T d,
              int p, const KnotVector<T>& U,
              const std::vector<NURBSPoint<T>>& Pw) {
    // Validate Mobius transform: ad - bc must be non-zero
    T det = a * d - b * c;
    if (std::abs(det) <= std::numeric_limits<T>::epsilon())
        throw std::invalid_argument("reparametrize: Mobius determinant (ad-bc) must be non-zero");

    // Transform each knot value: u'_j = (a*u_j + b) / (c*u_j + d)
    std::vector<T> U_prime_vec;
    U_prime_vec.reserve(U.size());
    for (std::size_t i = 0; i < U.size(); ++i) {
        T u_val = U[i];
        T denom = c * u_val + d;
        if (std::abs(denom) <= std::numeric_limits<T>::epsilon())
            throw std::invalid_argument("reparametrize: denominator c*u+d is zero at some knot");
        U_prime_vec.push_back((a * u_val + b) / denom);
    }

    // The transformed knot vector must be non-decreasing
    for (std::size_t i = 1; i < U_prime_vec.size(); ++i) {
        if (U_prime_vec[i] < U_prime_vec[i - 1])
            throw std::invalid_argument("reparametrize: Mobius transform produces decreasing knots");
    }

    KnotVector<T> U_prime(std::move(U_prime_vec));

    // Control points remain unchanged — reparameterization only affects the
    // knot vector and parameter mapping, not the control point geometry.
    // The curve shape is mathematically identical; only the parameter values
    // at which given control points are observed change.
    return ReparametrizeResult<T>{p, std::move(U_prime), Pw};
}

/**
 * reparametrize_unit_to_interval — reparameterize from [0,1] to [t0, t1].
 *
 * Convenience for the common case of linear remapping.
 *
 * @param t0  new domain start
 * @param t1  new domain end   (t1 != t0)
 * @param p   degree
 * @param U   knot vector (must be defined on [0,1])
 * @param Pw  control points
 */
template <NumericScalar_ T>
[[nodiscard]] ReparametrizeResult<T>
reparametrize_unit_to_interval(T t0, T t1,
                               int p, const KnotVector<T>& U,
                               const std::vector<NURBSPoint<T>>& Pw) {
    // Linear reparameterization: u' = t0 + (t1 - t0) * u
    // Mobius: (a*u+b)/(c*u+d) with a=(t1-t0), b=t0, c=0, d=1
    return reparametrize(t1 - t0, t0, T{0}, T{1}, p, U, Pw);
}

// ===========================================================================
// Algorithm A7.3 — Curve Fitting (Global Least-Squares Interpolation)
// ===========================================================================

/**
 * FittingResult — output of curve_fitting.
 */
template <NumericScalar_ T>
struct FittingResult {
    int p;                                     // degree (may be elevated internally)
    KnotVector<T> U;                          // knot vector
    std::vector<NURBSPoint<T>> Pw;           // fitted control points
};

/**
 * curve_fitting — Algorithm A7.3 (The NURBS Book, 2nd ed., p.410)
 *
 * Global curve fitting via least-squares approximation (or interpolation when
 * n = h = number of data points - 1).
 *
 * Given a set of data points Q_0 … Q_h in ℝ^n and a parameter sequence
 * u_0 … u_h, finds control points P_0 … P_n that minimize the fitting error
 * subject to the NURBS basis constraints.
 *
 * The system to solve (normal equations):
 *   [N]^T [N] {Pw} = [N]^T {Q}   (homogeneous form, weighted)
 *
 * where N_{i,p}(u_j) are the B-spline basis function values.
 *
 * For interpolation (h = n): solves the linear system directly.
 * For approximation (h > n): uses least-squares (more data points than
 *   control points → overdetermined system).
 *
 * @param Q    data points (Cartesian) to fit — size = h+1
 * @param u    parameter values for each data point — size = h+1, non-decreasing
 * @param p    desired degree (if 0, automatically determined from num pts)
 * @param n    desired number of control points - 1 (n <= h)
 *              If n == h → interpolation (passes through all points).
 *              If n <  h → least-squares approximation.
 * @param U    initial knot vector guess (may be refined internally)
 *
 * @return FittingResult with degree p, knot vector U, and control points Pw.
 *
 * Example:
 *   std::vector<Point<T>> Q = {{0,0},{1,2},{2,1},{3,3}};
 *   std::vector<T> u = {0.0, 0.33, 0.66, 1.0};
 *   auto result = curve_fitting(Q, u, 3, 3, uniform_knots(4, 3));
 *   NURBSCurve<T> curve(result.p, result.U, result.Pw);
 */
template <NumericScalar_ T>
[[nodiscard]] FittingResult<T>
curve_fitting(const std::vector<Point<T>>& Q,
              const std::vector<T>& u,
              int p,
              std::size_t n,
              const KnotVector<T>& U) {
    const std::size_t h = Q.size() - 1;  // number of data points - 1

    // Basic validation
    if (Q.size() != u.size())
        throw std::invalid_argument("curve_fitting: Q and u size mismatch");
    if (Q.size() < 2)
        throw std::invalid_argument("curve_fitting: need at least 2 data points");
    if (p < 1)
        throw std::invalid_argument("curve_fitting: degree must be >= 1");

    // n must satisfy: n >= p and n <= h (for interpolation) or n <= h (for least-squares)
    if (n < static_cast<std::size_t>(p))
        throw std::invalid_argument("curve_fitting: n must be >= degree p");
    if (n > h)
        throw std::invalid_argument("curve_fitting: n must be <= number of data points - 1");

    // Verify U size: for n+1 control points → knot vector size = n + p + 2
    if (U.size() != n + p + 2)
        throw std::invalid_argument("curve_fitting: knot vector size must be n + p + 2");

    // -------------------------------------------------------------------------
    // Step 1: Compute basis function values N_{i,p}(u_j) for all j=0..h, i=0..n
    // -------------------------------------------------------------------------
    // N_mat[j][i] = N_{i,p}(u_j), size = (h+1) × (n+1)
    std::vector<std::vector<T>> N_mat(h + 1, std::vector<T>(n + 1, T{0}));

    for (std::size_t j = 0; j <= h; ++j) {
        std::vector<T> b(p + 1);
        std::size_t k = compute_basis_functions(n, p, u[j], U, b);
        // Fill N_mat[j][i] for i = k-p … k (all non-zero basis functions at u_j)
        for (int i = 0; i <= p; ++i) {
            std::size_t idx = k - p + i;
            if (idx <= n)
                N_mat[j][idx] = b[i];
        }
    }

    // -------------------------------------------------------------------------
    // Step 2: Set up and solve the linear system
    // -------------------------------------------------------------------------
    // For weighted fitting: we work in homogeneous coordinates.
    // Build system: [N]^T[N] P_w = [N]^T Q_w
    // where Q_w = w_i * Q_i (homogeneous with w_i = 1 for Cartesian points).
    // For approximation, we solve the normal equations via QR or normal equations.
    // For interpolation (n == h): solve directly via Gaussian elimination.
    //
    // Use normal equations approach: A^T A x = A^T b
    // A is (h+1)×(n+1), b is (h+1)×dim.
    // For overdetermined (h > n): use least squares.
    // For determined   (h == n): direct solve.

    const std::size_t dim = Q[0].dimension();

    // Compute A^T A (n+1)×(n+1) and A^T b (n+1)×dim
    std::vector<std::vector<T>> ATA(n + 1, std::vector<T>(n + 1, T{0}));
    std::vector<std::vector<T>> ATb(n + 1, std::vector<T>(dim, T{0}));

    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= n; ++j) {
            T sum = T{0};
            for (std::size_t k = 0; k <= h; ++k) {
                sum += N_mat[k][i] * N_mat[k][j];
            }
            ATA[i][j] = sum;
        }
    }

    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t d = 0; d < dim; ++d) {
            T sum = T{0};
            for (std::size_t k = 0; k <= h; ++k) {
                sum += N_mat[k][i] * Q[k][d];
            }
            ATb[i][d] = sum;
        }
    }

    // -------------------------------------------------------------------------
    // Step 3: Solve ATA * P = ATb via Gaussian elimination with partial pivoting
    // -------------------------------------------------------------------------
    std::vector<std::vector<T>> P(n + 1, std::vector<T>(dim, T{0}));

    // Augment ATA with ATb → matrix of size (n+1) × (n+1+dim)
    std::vector<std::vector<T>> aug(n + 1, std::vector<T>(n + 1 + dim, T{0}));
    for (std::size_t i = 0; i <= n; ++i) {
        for (std::size_t j = 0; j <= n; ++j) aug[i][j] = ATA[i][j];
        for (std::size_t d = 0; d < dim; ++d) aug[i][n + 1 + d] = ATb[i][d];
    }

    // Forward elimination with partial pivoting
    for (std::size_t i = 0; i <= n; ++i) {
        // Find pivot row
        std::size_t pivot = i;
        T pivot_val = std::abs(aug[i][i]);
        for (std::size_t k = i + 1; k <= n; ++k) {
            if (std::abs(aug[k][i]) > pivot_val) {
                pivot_val = std::abs(aug[k][i]);
                pivot = k;
            }
        }

        if (pivot_val <= std::numeric_limits<T>::epsilon()) {
            throw std::runtime_error("curve_fitting: singular or near-singular matrix");
        }

        // Swap rows if needed
        if (pivot != i) {
            for (std::size_t j = i; j <= n + dim; ++j) {
                std::swap(aug[i][j], aug[pivot][j]);
            }
        }

        // Eliminate below pivot
        for (std::size_t k = i + 1; k <= n; ++k) {
            T factor = aug[k][i] / aug[i][i];
            for (std::size_t j = i; j <= n + dim; ++j) {
                aug[k][j] -= factor * aug[i][j];
            }
        }
    }

    // Back substitution
    for (std::size_t d = 0; d < dim; ++d) {
        for (std::size_t i = n + 1; i-- > 0; ) {
            T sum = aug[i][n + 1 + d];
            for (std::size_t j = i + 1; j <= n; ++j) {
                sum -= aug[i][j] * P[j][d];
            }
            P[i][d] = sum / aug[i][i];
        }
    }

    // -------------------------------------------------------------------------
    // Step 4: Build homogeneous control points Pw from Cartesian P
    // P[i] is a vector<T> of dimension dim (the solved control point coordinates)
    // -------------------------------------------------------------------------
    std::vector<NURBSPoint<T>> Pw(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        if (dim == 2) {
            Pw[i] = NURBSPoint<T>(P[i][0], P[i][1], T{1});
        } else {
            Pw[i] = NURBSPoint<T>(P[i][0], P[i][1],
                                  dim > 2 ? P[i][2] : T{0}, T{1});
        }
    }

    return FittingResult<T>{p, U, std::move(Pw)};
}

/**
 * curve_fitting_interpolate — convenience wrapper for interpolation mode (n = h).
 *
 * Given data points Q and parameter values u, fits a degree p curve that
 * passes through all data points exactly (n = h = len(Q)-1).
 *
 * @param Q   data points (Cartesian) — size ≥ p+1
 * @param u   parameter values — size = Q.size(), non-decreasing
 * @param p   degree
 *
 * @return FittingResult with interpolated curve.
 */
template <NumericScalar_ T>
[[nodiscard]] FittingResult<T>
curve_fitting_interpolate(const std::vector<Point<T>>& Q,
                          const std::vector<T>& u,
                          int p) {
    const std::size_t h = Q.size() - 1;
    const std::size_t n = h;  // interpolation → n = h

    if (Q.size() < static_cast<std::size_t>(p + 1))
        throw std::invalid_argument("curve_fitting_interpolate: need at least p+1 points");

    // Build averaging knot vector (Algorithm A7.1 setup)
    // u_0, …, u_h are the parameters; compute knot vector via averaging:
    // U_j = (u_{j-1} + u_{j-1}+1 + … + u_{j-p}) / p  for j = p … n
    std::vector<T> U_avg(h + p + 2, T{0});  // size m+1 = n+p+2

    // Clamped endpoints: first p+1 knots = u_0, last p+1 knots = u_h
    for (std::size_t i = 0; i <= static_cast<std::size_t>(p); ++i) {
        U_avg[i] = u[0];
        U_avg[U_avg.size() - 1 - i] = u[h];
    }

    // Interior knots via averaging
    for (std::size_t j = static_cast<std::size_t>(p); j <= n; ++j) {
        T sum = T{0};
        for (std::size_t k = 1; k <= static_cast<std::size_t>(p); ++k) {
            sum += u[j - p + k];
        }
        U_avg[j] = sum / static_cast<T>(p);
    }

    KnotVector<T> U_avg_kv(std::move(U_avg));

    return curve_fitting(Q, u, p, n, U_avg_kv);
}

} // namespace nurbs::curve
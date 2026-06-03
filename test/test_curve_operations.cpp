// test_curve_operations.cpp — Unit tests for Ch5 curve operations
// - curve_knot_insertion (Algorithm A5.4)
// - curve_degree_elevation (Algorithm A5.3)
// - curve_integrate (Algorithm A5.5)
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

#include "../include/nurbs/curve/curve_knot_insertion.hpp"
#include "../include/nurbs/curve/curve_degree_elevation.hpp"
#include "../include/nurbs/curve/curve_integrate.hpp"
#include "../include/nurbs/curve/nurbs_curve.hpp"

using namespace nurbs;
using namespace nurbs::core;
using namespace nurbs::curve;

// =============================================================================
// curve_knot_insertion — Algorithm A5.4
// =============================================================================

TEST_CASE("curve_knot_insertion: degree 1 line curve") {
    // Degree-1 B-spline with 3 control points: a straight line segment
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(2, 0, 1),
    };

    auto result = curve_knot_insertion(0.5, 1, U, Pw);

    CHECK(result.degree == 1);
    CHECK(result.knot_vector.size() == 5);     // one more knot
    CHECK(result.control_points.size() == 4);  // one more control point

    // Verify the new knot vector has 0.5 inserted
    CHECK(result.knot_vector[2] == 0.5);
}

TEST_CASE("curve_knot_insertion: preserves curve shape") {
    // Create a curve and verify that knot insertion doesn't change evaluated points
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);
    auto result = curve_knot_insertion(0.5, 2, U, Pw);

    NURBSCurve<double> crv2(result.degree, result.knot_vector, result.control_points);

    // Evaluate both curves at several parameters
    std::vector<double> params = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double u : params) {
        auto p1 = crv.evaluate(u);
        auto p2 = crv2.evaluate(u);
        CHECK(std::abs(p1.cart_x() - p2.cart_x()) < 1e-8);
        CHECK(std::abs(p1.cart_y() - p2.cart_y()) < 1e-8);
    }
}

TEST_CASE("curve_knot_insertion: knot already at max multiplicity") {
    // Inserting a knot that already exists at multiplicity p+1 should return unchanged
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};  // p=2 clamped
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto result = curve_knot_insertion(0.0, 2, U, Pw);  // 0 already at mult 3

    CHECK(result.degree == 2);
    CHECK(result.knot_vector.size() == U.size());  // unchanged
    CHECK(result.control_points.size() == Pw.size());
}

TEST_CASE("curve_insert_knot: convenience wrapper returns valid curve") {
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(2, 0, 1),
    };

    auto crv2 = curve_insert_knot(0.5, 1, U, Pw);

    CHECK(crv2.degree() == 1);
    CHECK(crv2.num_control_points() == 4);
    CHECK(crv2.num_knots() == 5);
}

// =============================================================================
// curve_degree_elevation — Algorithm A5.3
// =============================================================================

TEST_CASE("curve_degree_elevation: degree 1 to 2") {
    // Linear B-spline: degree 1 with 2 control points → elevate to degree 2
    // Requires knot vector {0,0,1,1} for degree 1 with 2 control points (n=1, m=n+p+1=1+1+1=3)
    // Actually: n=1, p=1, m=3 → 4 knots: {0,0,1,1} works for 2 CP (clamped)
    // But for degree elevation we need at least p+1=2 interior knots for a proper result.
    // Let's use degree-1 with 3 CPs: {0,0,0.5,1,1}
    std::vector<double> kv = {0, 0, 0.5, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 1, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto result = curve_degree_elevation(1, U, Pw, 1);

    CHECK(result.degree == 2);
    CHECK(result.knot_vector.size() > U.size());  // more knots after elevation
    CHECK(result.control_points.size() > Pw.size());
}

TEST_CASE("curve_degree_elevation: preserves curve shape") {
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 1, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);
    auto result = curve_degree_elevation(2, U, Pw, 1);

    NURBSCurve<double> crv_hi(result.degree, result.knot_vector, result.control_points);

    std::vector<double> params = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double u : params) {
        auto p1 = crv.evaluate(u);
        auto p2 = crv_hi.evaluate(u);
        CHECK(std::abs(p1.cart_x() - p2.cart_x()) < 1e-7);
        CHECK(std::abs(p1.cart_y() - p2.cart_y()) < 1e-7);
    }
}

TEST_CASE("curve_degree_elevation: invalid t throws") {
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 1, 1),
    };

    CHECK_THROWS(curve_degree_elevation(1, U, Pw, 0));   // t=0 invalid
    CHECK_THROWS(curve_degree_elevation(1, U, Pw, -1));  // t<0 invalid
}

TEST_CASE("curve_elevate_degree: convenience wrapper") {
    std::vector<double> kv = {0, 0, 0.5, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 1, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto crv2 = curve_elevate_degree(1, U, Pw, 1);

    CHECK(crv2.degree() == 2);
    CHECK(crv2.num_control_points() > 3);
}

// =============================================================================
// curve_integrate — Algorithm A5.5
// =============================================================================

TEST_CASE("curve_integrate: constant curve integrates to area") {
    // A degree-0 constant curve at y=1 from u=0 to u=1: S(u) = (1, 1, 0, 1)
    // Integral = ∫_0^1 (1,1,0,1) du = (1, 1, 0, 1)
    // For homogeneous NURBSPoint: w=1 so integral of (x,y,z,w) = (1,1,0,1)
    // Actually for B-spline degree 0, there's only one basis function = 1 over [0,1]
    // The integral of the curve point is just the point itself times the domain length.

    // Let's use a simple degree-1 line: C(u) = (u, 0, 0, 1) from 0 to 1
    // ∫_0^1 C(u) du = (0.5, 0, 0, 1)
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    // Control points: P0=(0,0), P1=(1,0) — line along x
    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto integral = curve_integrate(1, U, Pw);

    // For a straight line from (0,0) to (1,0), the area under is a triangle
    // with centroid at x=0.5. The integral should be (0.5, 0, 0, 1).
    CHECK(std::abs(integral.x() - 0.5) < 1e-8);
    CHECK(std::abs(integral.y()) < 1e-8);
    CHECK(std::abs(integral.w() - 1.0) < 1e-8);
}

TEST_CASE("curve_integrate: degree 2 Bezier segment") {
    // A degree-2 Bezier curve from u=0 to u=1 with control points
    // B0=(0,0), B1=(0.5,1), B2=(1,0). This is a parabola peaking at (0.5, 0.5).
    // For a quadratic Bezier, the integral is (1/3)*Σ B_i, scaled by interval length.
    // ∫_0^1 B(t)dt = (1/6) * (B0 + 4*B1 + B2) for standard parameterization
    // Actually: ∫_0^1 Σ C(2,i) t^i (1-t)^(2-i) dt = Σ C(2,i) * i!*(2-i)! / 3! = 1/3 * Σ C(2,i) / C(2,i)
    // = 1/3 * (B0 + B1 + B2) ???

    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 1, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto integral = curve_integrate(2, U, Pw);

    // For the quadratic Bezier, the exact integral in homogeneous form is:
    // (0.5, 1/3, 0, 1) — wait let me check the formula.
    // The Bernstein basis integral: ∫_0^1 B_{i,p}(t) dt = 1/(p+1)
    // So ∫_0^1 Σ P_i B_{i,p}(t) dt = Σ P_i * 1/(p+1)
    // For p=2: integral = (P0 + P1 + P2) / 3 = (0 + 0.5*1 + 1*0, 0+1*1+0, ...) / 3
    // = (0.5/3, 1/3, 0, 1) = (1/6, 1/3, 0, 1)
    CHECK(std::abs(integral.x() - 1.0/6.0) < 1e-8);
    CHECK(std::abs(integral.y() - 1.0/3.0) < 1e-8);
    CHECK(std::abs(integral.w() - 1.0) < 1e-8);
}

TEST_CASE("curve_integral_cartesian: simple line") {
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto integral_cart = curve_integral_cartesian(1, U, Pw);

    CHECK(std::abs(integral_cart[0] - 0.5) < 1e-8);
    CHECK(std::abs(integral_cart[1]) < 1e-8);
}

TEST_CASE("curve_integrate: multi-segment curve") {
    // A degree-2 curve with two Bezier segments: {0,0,0,0.5,0.5,1,1,1}
    std::vector<double> kv = {0, 0, 0, 0.5, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.25, 0.5, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(0.75, 0.5, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    auto integral = curve_integrate(2, U, Pw);

    // Sum of both segments: total integral = sum of individual segment integrals
    // Each segment has length 0.5 in parameter space
    // Segment 1: [0, 0.5], Segment 2: [0.5, 1]
    CHECK(integral.w() > 0.0);  // should have positive weight
    CHECK(std::abs(integral.x()) < 10.0);  // sanity bound
}
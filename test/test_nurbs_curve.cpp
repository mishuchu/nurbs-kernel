// test_nurbs_curve.cpp — Unit tests for NURBS curve (Ch5)
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

#include "../include/nurbs/curve/nurbs_curve.hpp"
#include "../include/nurbs/basis/knot_insertion.hpp"

using namespace nurbs;
using namespace nurbs::core;
using namespace nurbs::curve;

TEST_CASE("NURBSCurve: construction and basic accessors") {
    // Simple degree-2 curve with 4 control points (n=3)
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(0, 1, 1),
    };

    NURBSCurve<double> crv(2, std::move(U), std::move(Pw));

    CHECK(crv.degree() == 2);
    CHECK(crv.num_control_points() == 4);
    CHECK(crv.num_knots() == 7);
}

TEST_CASE("NURBSCurve: parameter_domain") {
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(0, 1, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);
    auto [u_min, u_max] = crv.parameter_domain();

    CHECK(u_min == 0.0);
    CHECK(u_max == 1.0);
}

TEST_CASE("NURBSCurve: endpoint evaluation") {
    // A degree-2 clamped B-spline should start at the first control point
    // and end at the last control point
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(0, 1, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);

    // At u=0, curve should start at P0 (with weight 1)
    auto C0 = crv.evaluate(0.0);
    CHECK(std::abs(C0.cart_x()) < 1e-8);
    CHECK(std::abs(C0.cart_y()) < 1e-8);

    // At u=1, curve should end at P3
    auto C1 = crv.evaluate(1.0);
    CHECK(std::abs(C1.cart_x()) < 1e-8);
    CHECK(std::abs(C1.cart_y() - 1.0) < 1e-8);
}

TEST_CASE("NURBSCurve: evaluate_derivatives returns correct size") {
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 0, 1),
        NURBSPoint<double>(1, 1, 1),
        NURBSPoint<double>(0, 1, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);
    auto ders = crv.evaluate_derivatives(0.5, 2);

    CHECK(ders.size() == 3); // C, C', C''
}

TEST_CASE("NURBSCurve: insert_knot preserves curve shape") {
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);

    // Insert knot at u=0.5
    auto crv2 = crv.insert_knot(0.5);

    CHECK(crv2.degree() == crv.degree());
    CHECK(crv2.num_control_points() == crv.num_control_points() + 1);

    // Evaluate at several points — should be numerically close
    std::vector<double> test_pts = {0.0, 0.25, 0.5, 0.75, 1.0};
    for (double u : test_pts) {
        auto p1 = crv.evaluate(u);
        auto p2 = crv2.evaluate(u);
        CHECK(std::abs(p1.cart_x() - p2.cart_x()) < 1e-8);
        CHECK(std::abs(p1.cart_y() - p2.cart_y()) < 1e-8);
    }
}

TEST_CASE("NURBSCurve: evaluate outside domain clamps") {
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    NURBSCurve<double> crv(2, U, Pw);

    // u < domain_min should clamp to P0
    auto Cneg = crv.evaluate(-0.5);
    auto C0 = crv.evaluate(0.0);
    CHECK(std::abs(Cneg.cart_x() - C0.cart_x()) < 1e-8);

    // u > domain_max should clamp to P2
    auto Chigh = crv.evaluate(2.0);
    auto C1 = crv.evaluate(1.0);
    CHECK(std::abs(Chigh.cart_x() - C1.cart_x()) < 1e-8);
}

TEST_CASE("NURBSCurve: construct_curve helper") {
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Pw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0.5, 1),
        NURBSPoint<double>(1, 1, 1),
    };

    auto crv = construct_curve(1, U, Pw);
    CHECK(crv.degree() == 1);
    CHECK(crv.num_control_points() == 3);
}
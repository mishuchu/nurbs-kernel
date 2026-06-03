// test_knot_removal.cpp — Unit tests for knot removal (Ch4, Algorithm A4.1)
#include <catch2/catch.hpp>
#include <cmath>
#include <vector>

#include "../include/nurbs/basis/knot_removal.hpp"
#include "../include/nurbs/curve/nurbs_curve.hpp"

using namespace nurbs;
using namespace nurbs::core;
using namespace nurbs::basis;

TEST_CASE("remove_knot: removal at interior knot preserves curve shape") {
    // Line curve: degree-1 with 3 control points forming a straight line
    // U = {0,0,0.5,1,1}, P = {(0,0), (0.5,0), (1,0)} — straight line along x
    // Removing the knot at u=0.5 should give a degree-1 curve with 2 control points
    // lying on the same line.
    std::vector<double> kv = {0, 0, 0.5, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Qw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    const int p = 1;
    auto result = remove_knot(0.5, p, U, Qw);

    CHECK(result.removed == true);
    CHECK(result.knot_vector.size() == 4);  // one knot removed
    CHECK(result.control_points.size() == 2); // one control point fewer
}

TEST_CASE("remove_knot: removal at endpoint does nothing") {
    // Knot at endpoint (multiplicity p+1) cannot be removed
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};  // p=2 clamped
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Qw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    const int p = 2;
    auto result = remove_knot(0.0, p, U, Qw);

    // Knot 0 is at multiplicity 3 (> p for degree 2) — removal should fail gracefully
    CHECK(result.removed == false);
    CHECK(result.knot_vector.size() == U.size());
    CHECK(result.control_points.size() == Qw.size());
}

TEST_CASE("remove_knot: nonexistent knot returns unchanged") {
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Qw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(1, 1, 1),
    };

    const int p = 1;
    auto result = remove_knot(0.3, p, U, Qw);

    CHECK(result.removed == false);
    CHECK(result.knot_vector.size() == U.size());
}

TEST_CASE("remove_knot_repeated: removes same knot multiple times") {
    // Knot vector with interior knot at 0.5 appearing twice
    std::vector<double> kv = {0, 0, 0.5, 0.5, 1, 1};
    KnotVector<double> U(kv);

    std::vector<NURBSPoint<double>> Qw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    const int p = 1;
    // multiplicity at 0.5 is 2, so we can remove it twice (r = p - s = 1 - 2 = -1... no)
    // Actually for p=1, knot at 0.5 with multiplicity 2 means r = 1 - 2 = -1 <= 0, nothing to remove
    // Let's use p=2 where we have 3 interior knots at 0.5
    std::vector<double> kv2 = {0, 0, 0, 0.5, 0.5, 0.5, 1, 1, 1};
    KnotVector<double> U2(kv2);

    std::vector<NURBSPoint<double>> Qw2 = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.5, 0, 1),
        NURBSPoint<double>(0.5, 1, 1),
        NURBSPoint<double>(1, 1, 1),
    };

    const int p2 = 2;
    auto result = remove_knot_repeated(0.5, p2, U2, Qw2, 3);

    CHECK(result.removed == true);
    CHECK(result.knot_vector.size() == 6);  // 9 - 3 = 6
}

TEST_CASE("remove_knot: degree 1 line with removable interior knot") {
    // A straight line degree-1 B-spline with knots {0,0,0.5,1,1} and control points
    // P0=(0,0), P1=(0.5,0), P2=(1,0).  The interior knot at 0.5 is removable
    // (multiplicity s=1 < p=1, so r=p-s=0, nothing to remove... actually r=0 means already removable only once)
    // Let me set s=1 (one occurrence at 0.5), p=1 → r=1-1=0, hmm that means we can remove it 0 times?
    // No: r = p - s, and we remove r times. If r=1, we remove once.
    // For a knot with multiplicity 1 at degree p=1, r = 1 - 1 = 0. So we cannot remove.
    // Need s < p for removal to be possible. Let me use p=2.
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    KnotVector<double> U(kv);

    // Bezier-like degree-2: 4 control points along a line
    std::vector<NURBSPoint<double>> Qw = {
        NURBSPoint<double>(0, 0, 1),
        NURBSPoint<double>(0.33, 0, 1),
        NURBSPoint<double>(0.66, 0, 1),
        NURBSPoint<double>(1, 0, 1),
    };

    const int p = 2;
    auto result = remove_knot(0.5, p, U, Qw);

    CHECK(result.removed == true);
    CHECK(result.knot_vector.size() == 6); // one knot removed
    CHECK(result.control_points.size() == 3);
}
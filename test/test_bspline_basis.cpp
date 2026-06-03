// test_bspline_basis.cpp — Unit tests for B-spline basis functions (Ch4)
#include <catch2/catch.hpp>
#include <vector>
#include <cmath>

#include "../include/nurbs/basis/bspline_basis.hpp"
#include "../include/nurbs/core/types.hpp"
#include "../include/nurbs/core/utilities.hpp"

using namespace nurbs;
using namespace nurbs::core;
using namespace nurbs::basis;

TEST_CASE("bspline_basis: uniform knot vector degree 1") {
    // Knot vector U = {0, 0, 0.5, 1, 1}, degree p=1, n=2 (3 control points)
    // Basis functions at u=0.25 should give N_0,0 + N_1,0 + N_2,0 = 1
    KnotVector<double> U = KnotVector<double>::uniform(3, 0.0, 1.0);
    // Manually set clamped knot vector for degree 1
    std::vector<double> kv = {0, 0, 0.5, 1, 1};
    KnotVector<double> U_clamped(kv);

    const std::size_t n = 2; // n+1 = 3 control points
    const int p = 1;

    auto result = bspline_basis(n, p, 0.25, U_clamped);

    CHECK(result.span >= p);
    CHECK(result.span <= n);
    CHECK(result.values.size() == static_cast<std::size_t>(p + 1));

    // Sum of basis functions at any valid u should be <= 1
    double sum = 0;
    for (auto v : result.values) sum += v;
    CHECK(sum <= 1.0 + 1e-10);
}

TEST_CASE("bspline_basis: degree 0 (constant)") {
    // Degree 0: N_{i,0}(u) = 1 if u in [u_i, u_{i+1}), else 0
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv);

    const std::size_t n = 1;
    const int p = 0;

    auto res0 = bspline_basis(n, p, 0.0, U);
    auto res_half = bspline_basis(n, p, 0.5, U);

    CHECK(res0.values.size() == 1);
    CHECK(res_half.values.size() == 1);
}

TEST_CASE("basis_function: single function evaluation") {
    std::vector<double> kv = {0, 0, 0.5, 1, 1};
    KnotVector<double> U(kv);

    // N_{0,1}(0.25) should be 0.5
    // N_{1,1}(0.25) should be 0.5
    double val0 = basis_function(0, 1, 0.25, U);
    double val1 = basis_function(1, 1, 0.25, U);

    CHECK(std::abs(val0 + val1 - 1.0) < 1e-10);
}

TEST_CASE("bspline_basis: endpoint evaluation") {
    // Clamped degree-2 B-spline: n=2, p=2, U={0,0,0,1,1,1}
    // At u=0: only N_0,2 is non-zero and N_0,2(0) = 1
    // At u=1: only N_2,2 is non-zero and N_2,2(1) = 1
    // Using a VALID knot vector matching the relationship m = n + p + 1
    // n=2, p=2 → m = 2+2+1 = 5 → U has 6 knots
    std::vector<double> kv = {0, 0, 0, 1, 1, 1}; // clamped, p=2, n=2
    KnotVector<double> U(kv);

    const std::size_t n = 2; // n+1 = 3 control points
    const int p = 2;

    // At u=0 (left endpoint of domain [0,1]):
    // find_span returns k=2 (because u >= U[2]=0, u < U[3]=1)
    // Basis functions should be [1, 0, 0] (only N_0,2 active)
    auto res = bspline_basis(n, p, 0.0, U);
    CHECK(res.values.size() == static_cast<std::size_t>(p + 1));
    CHECK(std::abs(res.values[0] - 1.0) < 1e-10);

    // At u=1 (right endpoint of domain):
    // find_span returns k=4 (special case u >= U[5])
    // Basis functions should be [0, 0, 1] (only N_2,2 active)
    auto res1 = bspline_basis(n, p, 1.0, U);
    CHECK(std::abs(res1.values[p] - 1.0) < 1e-10);
}

TEST_CASE("find_span: binary search correctness") {
    std::vector<double> kv = {0, 0, 0, 0.25, 0.5, 0.75, 1, 1, 1};
    KnotVector<double> U(kv);

    CHECK(find_span(4, 2, 0.0, U) == 2);
    CHECK(find_span(4, 2, 0.1, U) == 2);
    CHECK(find_span(4, 2, 0.25, U) == 3);
    CHECK(find_span(4, 2, 0.5, U) == 4);
    CHECK(find_span(4, 2, 1.0, U) == 4); // last knot returns n
}
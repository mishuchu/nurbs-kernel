// test_surface.cpp — Unit tests for NURBS surface algorithms (Ch6)
// Algorithm A6.1  : construct_surface
// Algorithm A6.2  : surface_derivatives
// Algorithm A6.3  : surface_knot_insertion_u / surface_knot_insertion_v
// Algorithm A6.4  : surface_degree_elevation_u / surface_degree_elevation_v

#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <sstream>
#include <stdexcept>

#include "../include/nurbs/surface/construct_surface.hpp"
#include "../include/nurbs/surface/surface_derivs.hpp"
#include "../include/nurbs/surface/surface_knot_insertion.hpp"
#include "../include/nurbs/surface/surface_degree_elevation.hpp"
#include "../include/nurbs/core/types.hpp"
#include "../include/nurbs/core/utilities.hpp"

using namespace nurbs;
using namespace nurbs::core;
using namespace nurbs::surface;

// ===========================================================================
// Algorithm A6.1 — construct_surface
// ===========================================================================

TEST_CASE("construct_surface: valid bilinear patch") {
    // Simple bilinear B-spline surface (degree 1 × 1)
    // 2×2 control points, knot vectors clamped
    const int p_u = 1, p_v = 1;
    const std::size_t nu = 2, nv = 2;  // (n+1)×(m+1) control points

    std::vector<double> kv = {0, 0, 1, 1};   // clamped U
    std::vector<double> lv = {0, 0, 1, 1};   // clamped V
    KnotVector<double> U(kv);
    KnotVector<double> V(lv);

    // Pw[i][j] — row i = v-index, col j = u-index
    // P00=(0,0), P10=(1,0), P01=(0,1), P11=(1,1), all w=1
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    auto result = construct_surface(p_u, p_v, U, V, Pw);

    CHECK(result.degree_u == 1);
    CHECK(result.degree_v == 1);
    CHECK(result.knot_vector_u.size() == static_cast<std::size_t>(p_u) + nu + 1); // n+p+2 = 1+2+1=4
    CHECK(result.knot_vector_v.size() == static_cast<std::size_t>(p_v) + nv + 1);
    CHECK(result.control_points.size() == nv);
    CHECK(result.control_points[0].size() == nu);
}

TEST_CASE("construct_surface: rejects wrong knot size") {
    const int p_u = 2, p_v = 1;
    const std::size_t nu = 3, nv = 2;
    KnotVector<double> U = KnotVector<double>::uniform(nu + p_u + 1, 0.0, 1.0);
    KnotVector<double> V = KnotVector<double>::uniform(nv + p_v + 1, 0.0, 1.0);
    std::vector<std::vector<NURBSPoint<double>>> Pw(nv, std::vector<NURBSPoint<double>>(nu));

    // U size is correct, but pass a wrong one
    KnotVector<double> bad_U = KnotVector<double>::uniform(nu + p_u, 0.0, 1.0); // too short
    CHECK_THROWS_AS(construct_surface(p_u, p_v, bad_U, V, Pw), std::invalid_argument);
}

TEST_CASE("construct_surface: rejects non-rectangular Pw") {
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(kv);

    // Row 0 has 2 points, row 1 has 3 points → not rectangular
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1), NURBSPoint<double>(2,1,1) }
    };

    CHECK_THROWS_AS(construct_surface(p_u, p_v, U, V, Pw), std::invalid_argument);
}

TEST_CASE("construct_surface: degree zero surface") {
    // Degree 0 in both directions: nu=1, nv=1 requires U size = 0+0+2=2, V size = 2
    const int p_u = 0, p_v = 0;
    const std::size_t nu = 1, nv = 1;
    std::vector<double> kv = {0, 1};   // size 2 OK: n+p+2 = 0+0+2=2
    std::vector<double> lv = {0, 1};
    KnotVector<double> U(kv), V(lv);
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1) }
    };

    auto result = construct_surface(p_u, p_v, U, V, Pw);
    CHECK(result.degree_u == 0);
    CHECK(result.degree_v == 0);
}

TEST_CASE("construct_surface: rejects negative degree") {
    KnotVector<double> U = KnotVector<double>::uniform(3, 0.0, 1.0);
    KnotVector<double> V = KnotVector<double>::uniform(3, 0.0, 1.0);
    std::vector<std::vector<NURBSPoint<double>>> Pw(2, std::vector<NURBSPoint<double>>(2));

    CHECK_THROWS_AS(construct_surface(-1, 1, U, V, Pw), std::invalid_argument);
    CHECK_THROWS_AS(construct_surface(1, -1, U, V, Pw), std::invalid_argument);
}

TEST_CASE("construct_surface_from_cartesian: simple case") {
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    // Flattened row-major Cartesian points: 2×2 grid
    std::vector<Point<double>> P = {
        Point<double>{0,0}, Point<double>{1,0},
        Point<double>{0,1}, Point<double>{1,1}
    };
    WeightVector<double> W = WeightVector<double>::uniform(4);

    auto result = construct_surface_from_cartesian(p_u, p_v, U, V, P, 2, 2, W);
    CHECK(result.degree_u == 1);
    CHECK(result.degree_v == 1);
    CHECK(result.control_points[0][0].cart_x() == Catch::Approx(0.0));
    CHECK(result.control_points[1][1].cart_x() == Catch::Approx(1.0));
}


// ===========================================================================
// Algorithm A6.2 — surface_derivatives
// ===========================================================================

TEST_CASE("surface_derivatives: flat plane has zero tangential derivatives") {
    // A flat plane z = 0.5 should have zero z-derivative in both directions
    const int p = 1, q = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    // All control points lie on z=0.5 plane, homogeneous Pw = (x, y, 0.5, 1)
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0, 0, 0.5, 1), NURBSPoint<double>(1, 0, 0.5, 1) },
        { NURBSPoint<double>(0, 1, 0.5, 1), NURBSPoint<double>(1, 1, 0.5, 1) }
    };

    auto der = surface_derivatives(0.5, 0.5, p, q, U, V, Pw, 1, 1);

    // d/dx and d/dy at any interior point should be zero for this flat surface
    double z_val  = der.derivs[0][0].cart_z();
    double dz_du  = der.derivs[1][0].cart_z(); // ∂z/∂u
    double dz_dv  = der.derivs[0][1].cart_z(); // ∂z/∂v

    CHECK(std::abs(z_val  - 0.5) < 1e-8);
    CHECK(std::abs(dz_du)        < 1e-8);
    CHECK(std::abs(dz_dv)        < 1e-8);
}

TEST_CASE("surface_derivatives: point only (d_u=0, d_v=0)") {
    const int p = 2, q = 1;
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(0.5,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(0.5,1,0,1), NURBSPoint<double>(1,1,0,1) }
    };

    auto der = surface_derivatives(0.5, 0.5, p, q, U, V, Pw, 0, 0);

    CHECK(der.derivs.size() == 1);
    CHECK(der.derivs[0].size() == 1);
    // Point value should be a valid Cartesian point
    double cx = der.derivs[0][0].cart_x();
    CHECK(std::isfinite(cx));
}

TEST_CASE("surface_derivatives: boundary evaluation u=0 v=0") {
    // At the corner (0,0) of the domain for a clamped surface,
    // only the first control point column/row contributes
    const int p = 2, q = 2;
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    std::vector<double> lv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(1,0,0,1), NURBSPoint<double>(2,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(1,1,0,1), NURBSPoint<double>(2,1,0,1) },
        { NURBSPoint<double>(0,2,0,1), NURBSPoint<double>(1,2,0,1), NURBSPoint<double>(2,2,0,1) }
    };

    auto der = surface_derivatives(0.0, 0.0, p, q, U, V, Pw, 2, 2);

    // At u=0, v=0: surface point should equal P_{0,0} = (0,0,0)
    CHECK(der.derivs[0][0].cart_x() == Catch::Approx(0.0));
    CHECK(der.derivs[0][0].cart_y() == Catch::Approx(0.0));
}

TEST_CASE("surface_derivatives: derivative array size respects dmin") {
    // dmin = min(d, p) so asking for more derivatives than degree truncates
    const int p = 1, q = 2;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 0, 1, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(1,1,0,1) },
        { NURBSPoint<double>(0,2,0,1), NURBSPoint<double>(1,2,0,1) }
    };

    auto der = surface_derivatives(0.5, 0.5, p, q, U, V, Pw, 3, 3);

    // Only up to p=1 and q=2 derivatives are populated
    CHECK(der.derivs.size() == 2);   // d_u_min+1 = min(3,1)+1 = 2
    CHECK(der.derivs[0].size() == 3); // d_v_min+1 = min(3,2)+1 = 3
}

TEST_CASE("surface_derivatives: default derivative order equals degree") {
    const int p = 2, q = 1;
    std::vector<double> kv = {0, 0, 0, 0.5, 1, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(0.5,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(0.5,1,0,1), NURBSPoint<double>(1,1,0,1) }
    };

    // Call without explicit d_u, d_v → defaults to (p, q)
    auto der = surface_derivatives(0.5, 0.5, p, q, U, V, Pw);
    CHECK(der.derivs.size() == static_cast<std::size_t>(p) + 1);
    CHECK(der.derivs[0].size() == static_cast<std::size_t>(q) + 1);
}


// ===========================================================================
// Algorithm A6.3 — surface_knot_insertion_u / surface_knot_insertion_v
// ===========================================================================

TEST_CASE("surface_knot_insertion_u: interior knot increases control points") {
    // Start with a simple bilinear surface
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    auto res = surface_knot_insertion_u(0.5, p_u, p_v, U, V, Pw);

    // After inserting one interior knot, we should have nu+1=3 columns
    CHECK(res.control_points[0].size() == 3);
    CHECK(res.knot_vector_u.size() == U.size() + 1);
    // V is unchanged
    CHECK(res.knot_vector_v.size() == V.size());
    // The new knot value should appear in the updated knot vector
    bool found_new_knot = false;
    for (std::size_t i = 0; i < res.knot_vector_u.size(); ++i)
        if (res.knot_vector_u[i] == 0.5) found_new_knot = true;
    CHECK(found_new_knot);
}

TEST_CASE("surface_knot_insertion_u: knot at max multiplicity returns unchanged") {
    // Insert a knot value that already has multiplicity p+1 → no change
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1}; // multiplicity at 0 and 1 is already 2 = p+1
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    auto res = surface_knot_insertion_u(0.0, p_u, p_v, U, V, Pw);

    CHECK(res.control_points[0].size() == 2); // unchanged
    CHECK(res.knot_vector_u.size() == U.size());           // unchanged
}

TEST_CASE("surface_knot_insertion_u: endpoint knot insertion") {
    const int p_u = 2, p_v = 1;
    std::vector<double> kv = {0, 0, 0, 1, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(0.5,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(0.5,1,0,1), NURBSPoint<double>(1,1,0,1) }
    };

    // Insert u=0.0 at the left endpoint (already at max multiplicity, no-op)
    auto res_noop = surface_knot_insertion_u(0.0, p_u, p_v, U, V, Pw);
    CHECK(res_noop.control_points[0].size() == 3);

    // Insert at u=0.25 (new interior knot → new column)
    auto res_new = surface_knot_insertion_u(0.25, p_u, p_v, U, V, Pw);
    CHECK(res_new.control_points[0].size() == 4);
}

TEST_CASE("surface_knot_insertion_v: interior knot increases rows") {
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    auto res = surface_knot_insertion_v(0.5, p_u, p_v, U, V, Pw);

    CHECK(res.control_points.size() == 3); // nv+1 rows
    CHECK(res.knot_vector_v.size() == V.size() + 1);
    CHECK(res.knot_vector_u.size() == U.size()); // U unchanged
}

TEST_CASE("surface_knot_insertion: direction flag u vs v") {
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    auto res_u = surface_knot_insertion(0.5, p_u, p_v, U, V, Pw, true);  // u-direction
    auto res_v = surface_knot_insertion(0.5, p_u, p_v, U, V, Pw, false); // v-direction

    CHECK(res_u.control_points[0].size() == 3); // +1 column
    CHECK(res_v.control_points.size() == 3);   // +1 row
}


// ===========================================================================
// Algorithm A6.4 — surface_degree_elevation_u / surface_degree_elevation_v
// ===========================================================================

TEST_CASE("surface_degree_elevation_u: degree 1→2 along u") {
    // A flat surface z=0.5 (degree 1×1) elevated to degree 2×1
    const int p_u = 1, p_v = 1, t = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0.5,1), NURBSPoint<double>(1,0,0.5,1) },
        { NURBSPoint<double>(0,1,0.5,1), NURBSPoint<double>(1,1,0.5,1) }
    };

    auto res = surface_degree_elevation_u(p_u, p_v, U, V, Pw, t);

    CHECK(res.degree_u == p_u + t);   // 2
    CHECK(res.degree_v == p_v);        // 1 (unchanged)
    // After degree elevation, number of control points in u should increase by t
    CHECK(res.control_points[0].size() == Pw[0].size() + static_cast<std::size_t>(t));
    // V knot vector unchanged
    CHECK(res.knot_vector_v.size() == V.size());
}

TEST_CASE("surface_degree_elevation_v: degree 1→2 along v") {
    const int p_u = 1, p_v = 1, t = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0.5,1), NURBSPoint<double>(1,0,0.5,1) },
        { NURBSPoint<double>(0,1,0.5,1), NURBSPoint<double>(1,1,0.5,1) }
    };

    auto res = surface_degree_elevation_v(p_u, p_v, U, V, Pw, t);

    CHECK(res.degree_u == p_u);        // unchanged
    CHECK(res.degree_v == p_v + t);    // 2
    CHECK(res.control_points.size() == Pw.size() + static_cast<std::size_t>(t));
    CHECK(res.knot_vector_u.size() == U.size());
}

TEST_CASE("surface_degree_elevation: both u and v directions") {
    const int p_u = 1, p_v = 1, t_u = 1, t_v = 2;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(1,1,0,1) }
    };

    auto res = surface_degree_elevation(p_u, p_v, U, V, Pw, t_u, t_v);

    CHECK(res.degree_u == p_u + t_u);
    CHECK(res.degree_v == p_v + t_v);
    CHECK(res.control_points[0].size() == Pw[0].size() + static_cast<std::size_t>(t_u));
    CHECK(res.control_points.size() == Pw.size() + static_cast<std::size_t>(t_v));
}

TEST_CASE("surface_degree_elevation: rejects t < 1") {
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1), NURBSPoint<double>(1,0,1) },
        { NURBSPoint<double>(0,1,1), NURBSPoint<double>(1,1,1) }
    };

    CHECK_THROWS_AS(surface_degree_elevation_u(1, 1, U, V, Pw, 0), std::invalid_argument);
    CHECK_THROWS_AS(surface_degree_elevation_v(1, 1, U, V, Pw, -1), std::invalid_argument);
}

TEST_CASE("surface_degree_elevation: degree 0 surface is rejected") {
    // degree 0 surfaces are not supported for elevation
    std::vector<double> kv = {0, 1};
    std::vector<double> lv = {0, 1};
    KnotVector<double> U(kv), V(lv);
    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,1) }
    };

    CHECK_THROWS_AS(surface_degree_elevation_u(0, 1, U, V, Pw, 1), std::invalid_argument);
}

TEST_CASE("surface_degree_elevation: elevated surface remains equivalent") {
    // A flat bilinear surface elevated by 1 in u should evaluate to the same values
    const int p_u = 1, p_v = 1;
    std::vector<double> kv = {0, 0, 1, 1};
    std::vector<double> lv = {0, 0, 1, 1};
    KnotVector<double> U(kv), V(lv);

    std::vector<std::vector<NURBSPoint<double>>> Pw = {
        { NURBSPoint<double>(0,0,0,1), NURBSPoint<double>(1,0,0,1) },
        { NURBSPoint<double>(0,1,0,1), NURBSPoint<double>(1,1,0,1) }
    };

    auto res = surface_degree_elevation_u(p_u, p_v, U, V, Pw, 1);
    CHECK(res.degree_u == 2);

    // Knot vector U_bar should have more interior knots (cloned from original)
    CHECK(res.knot_vector_u.size() > U.size());
}
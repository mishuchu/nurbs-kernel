// nurbs.hpp — Master include for the NURBS kernel library
#pragma once

// Core
#include "nurbs/core/concepts.hpp"
#include "nurbs/core/types.hpp"
#include "nurbs/core/numeric.hpp"
#include "nurbs/core/utilities.hpp"

// Basis functions (Ch4)
#include "nurbs/basis/bspline_basis.hpp"
#include "nurbs/basis/knot_insertion.hpp"
#include "nurbs/basis/knot_refinement.hpp"
#include "nurbs/basis/degree_elevation.hpp"

// Curves (Ch5)
#include "nurbs/curve/nurbs_curve.hpp"
#include "nurbs/curve/curve_derivatives.hpp"
#include "nurbs/curve/curve_inversion.hpp"

// Surfaces (Ch6) — structural stubs
#include "nurbs/surface/nurbs_surface.hpp"
// nurbs_surface.hpp — NURBS surface class (The NURBS Book, Ch6)
// Placeholder stub — full implementation is Ch6 milestone
#pragma once

#include "../core/concepts.hpp"
#include "../core/types.hpp"
#include "../core/utilities.hpp"

namespace nurbs::surface {

// Forward-declare surface concept
template <nurbs::core::NumericScalar_ T = double>
class NURBSSurface {
public:
    using scalar_type = T;

    NURBSSurface() = default;
    // TODO: full Ch6 implementation
};

} // namespace nurbs::surface
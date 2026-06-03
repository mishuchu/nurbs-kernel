# Catch2 amalgamated (header-only) — no find_package needed
set(CATCH2_AMALGAMATED_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../3rdparty/catch2")

# ------------------------------------------------------------------ test executable
add_executable(tests
    test/test_main.cpp
    test/test_bspline_basis.cpp
    test/test_nurbs_curve.cpp
)
target_link_libraries(tests PRIVATE nurbs-kernel)
target_compile_options(tests PRIVATE -Wall -Wextra -std=c++20)

# Register with CTest
include(CTest)
register_tests(TARGET tests)
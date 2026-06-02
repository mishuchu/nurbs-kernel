cmake_minimum_required(VERSION 3.20)
project(nurbs-kernel VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} PRIVATE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(${PROJECT_NAME} PRIVATE /W4)
endif()

# Header-only: no target_sources needed, just PUBLIC header dirs
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)

# Alias for modern usage
add_library(${PROJECT_NAME}::nurbs-kernel ALIAS ${PROJECT_NAME})

# ----- Testing -----
include(CTest)
if(BUILD_TESTING)
    find_package(Catch2 QUIET)
    if(Catch2_FOUND)
        add_executable(tests
            test/test_main.cpp
            test/test_bspline_basis.cpp
            test/test_nurbs_curve.cpp
        )
        target_link_libraries(tests PRIVATE Catch2::Catch2 nurbs-kernel)
        target_compile_options(tests PRIVATE -Wall -Wextra)
        register_tests(TARGET tests)
    endif()
endif()

# ----- Installation -----
include(GNUInstallDirs)
install(DIRECTORY include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}Targets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(EXPORT ${PROJECT_NAME}Targets
    FILE ${PROJECT_NAME}Targets.cmake
    NAMESPACE ${PROJECT_NAME}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
)

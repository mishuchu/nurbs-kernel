cmake_minimum_required(VERSION 3.20)
project(nurbs-kernel VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# ----- CTest must be included BEFORE any targets -----
include(CTest)

# ----- Library: header-only INTERFACE -----
add_library(${PROJECT_NAME} INTERFACE)
target_include_directories(${PROJECT_NAME} INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(${PROJECT_NAME} INTERFACE cxx_std_20)
add_library(${PROJECT_NAME}::nurbs-kernel ALIAS ${PROJECT_NAME})

# Warnings
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(${PROJECT_NAME} INTERFACE -Wall -Wextra -Wpedantic)
elseif(MSVC)
    target_compile_options(${PROJECT_NAME} INTERFACE /W4)
endif()

# ----- Testing: amalgamated Catch2 -----
if(BUILD_TESTING)
    # Fetch Catch2 amalgamated if not present
    set(CATCH2_AMALGAMATED_DIR "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/catch2")
    if(NOT EXISTS "${CATCH2_AMALGAMATED_DIR}/catch_amalgamated.hpp")
        message(STATUS "Fetching Catch2 amalgamated from GitHub...")
        file(MAKE_DIRECTORY "${CATCH2_AMALGAMATED_DIR}")
        file(DOWNLOAD
            "https://raw.githubusercontent.com/catchorg/Catch2/devel/extras/catch_amalgamated.cpp"
            "${CATCH2_AMALGAMATED_DIR}/catch_amalgamated.cpp"
            SHOW_PROGRESS
        )
        file(DOWNLOAD
            "https://raw.githubusercontent.com/catchorg/Catch2/devel/extras/catch_amalgamated.hpp"
            "${CATCH2_AMALGAMATED_DIR}/catch_amalgamated.hpp"
            SHOW_PROGRESS
        )
    endif()

    if(NOT EXISTS "${CATCH2_AMALGAMATED_DIR}/catch_amalgamated.hpp")
        message(FATAL_ERROR "3rdparty/catch2/catch_amalgamated.hpp not found")
    endif()
    message(STATUS "Catch2: using amalgamated headers")

    # Compile the amalgamated .cpp as an OBJECT library
    add_library(Catch2_amalgamated OBJECT
        "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/catch2/catch_amalgamated.cpp"
    )
    target_include_directories(Catch2_amalgamated PUBLIC
        "${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/catch2"
    )
    target_compile_options(Catch2_amalgamated PRIVATE -Wall -Wextra -std=c++20)
    target_compile_definitions(Catch2_amalgamated PUBLIC CATCH_CONFIG_MAIN)
    add_library(Catch2::Catch2 ALIAS Catch2_amalgamated)

    # Build the test executable
    add_executable(tests
        test/test_main.cpp
        test/test_bspline_basis.cpp
        test/test_nurbs_curve.cpp
        test/test_surface.cpp
    )
    target_link_libraries(tests PRIVATE Catch2::Catch2 nurbs-kernel)
    target_compile_options(tests PRIVATE -Wall -Wextra -std=c++20)
    add_test(NAME tests COMMAND tests)
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
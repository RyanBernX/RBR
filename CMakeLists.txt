cmake_minimum_required(VERSION 3.5)

project(RBR VERSION 1.0.0)

# disable in-source builds
set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

# set C99 standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set C++11 standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set module path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

# define options
option(BUILD_SHARED_LIBS "Enable shared library. Default: ON" ON)
option(BUILD_MATLAB_INTERFACE "Enable the build of matlab interface. Default: OFF" OFF)
option(USE_ILP64 "Build RBR with 64-bit integer interface. Default: OFF" OFF)

# build types
if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE is unset, defaulting to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# compilers support
if ("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    add_compile_options(-Wall -Wextra -Wpedantic)
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        add_compile_options(-g)
    endif()
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    # note: not working for current version because C++17 support is poor
    add_compile_options(-Wall -Wextra -inline-forceinline)
elseif ("${CMAKE_C_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_C_COMPILER_ID}" STREQUAL "AppleClang")
    add_compile_options(-Wall -Wextra -Wpedantic)
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        add_compile_options(-g)
    elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
        # C99 是允许使用的，但 Clang 依然可能会给出警告，因此屏蔽
        add_compile_options(-Wno-c99-extensions)
    endif()
else()
    message(WARNING "Unsupported C compiler ${CMAKE_C_COMPILER_ID}")
endif()

# MacOS extra search paths
if (APPLE AND ENABLE_EXTRA_PREFIX_PATH_MACOS)
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};/opt/homebrew/opt/openblas;/opt/homebrew/opt/libomp")
    # try to find frameworks after standard libraries or headers
    set(CMAKE_FIND_FRAMEWORK "LAST")
endif()

# openmp dep
find_package(OpenMP REQUIRED COMPONENTS C)

# blas dep
find_package(CBLAS REQUIRED)

# check whether it is MKL
if (BLAS_LIBRARIES MATCHES "mkl")
    message(STATUS "Enable MKL support in RBR.")
    add_definitions(-DHAVE_MKL)
endif()

# ilp64
if (USE_ILP64)
    add_compile_options(-DRBR_ILP64)
endif()

# subdirs
add_subdirectory(src)

# output compile summary
include(cmake/utils.cmake)
message("-- Configuration summary for RBR:")
message("   -- prefix: ${CMAKE_INSTALL_PREFIX}")
message("   -- build: ${CMAKE_BUILD_TYPE}")
message("   -- SHARED_LIBS: ${BUILD_SHARED_LIBS}")
message("   -- ILP64: ${USE_ILP64}")
get_target_property(MAIN_FLAGS RBR COMPILE_OPTIONS)
string(REPLACE ";" " " MAIN_FLAGS_STR "${MAIN_FLAGS}")
cprsummary("CC" "${CMAKE_C_COMPILER}"
    "${CMAKE_C_FLAGS_DEBUG}"
    "${CMAKE_C_FLAGS_MINSIZEREL}"
    "${CMAKE_C_FLAGS_RELEASE}"
    "${CMAKE_C_FLAGS_RELWITHDEBINFO}"
    "${MAIN_FLAGS_STR}")

cprsummary("CXX" "${CMAKE_CXX_COMPILER}"
    "${CMAKE_CXX_FLAGS_DEBUG}"
    "${CMAKE_CXX_FLAGS_MINSIZEREL}"
    "${CMAKE_CXX_FLAGS_RELEASE}"
    "${CMAKE_CXX_FLAGS_RELWITHDEBINFO}"
    "${MAIN_FLAGS_STR}")

libsummary("CBLAS" "${CBLAS_INCLUDE_DIRS}" "${CBLAS_LIBRARIES}")

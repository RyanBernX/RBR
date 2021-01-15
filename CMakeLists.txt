cmake_minimum_required(VERSION 3.3)

project(RBR VERSION 0.1.0)

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

file(GLOB_RECURSE C_SOURCE ${PROJECT_SOURCE_DIR}/src/*.c)
file(GLOB_RECURSE CXX_SOURCE ${PROJECT_SOURCE_DIR}/src/*.cpp)
list(FILTER C_SOURCE EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/main\.c")
list(FILTER C_SOURCE EXCLUDE REGEX "${PROJECT_SOURCE_DIR}/src/mex_rbr\.c")
set(SOURCE ${C_SOURCE} ${CXX_SOURCE})
set(MATLAB_SOURCE ${PROJECT_SOURCE_DIR}/src/mex_rbr.c)

# define options
option(BUILD_SHARED_LIBS "Enable shared library. Default: ON" ON)
option(BUILD_MATLAB_INTERFACE
    "Enable the build of matlab interface. Default: OFF" OFF)

# define targets
add_executable(rbr src/main.c)
add_library(RBR ${SOURCE})
target_link_libraries(rbr RBR)

configure_file(include/config.h.in config.h)

if ("${CMAKE_BUILD_TYPE}" STREQUAL "")
    message(STATUS "CMAKE_BUILD_TYPE is unset, defaulting to Debug")
    set(CMAKE_BUILD_TYPE "Debug")
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -m64 -march=native")
    if ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
        add_compile_options(-Wpedantic)
    endif()
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    # note: not working for current version because C++17 support is poor
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -inline-forceinline")
else()
    message(WARNING "Unsupport CXX compiler ${CMAKE_CXX_COMPILER_ID}")
endif()

include_directories("${PROJECT_BINARY_DIR}" "${PROJECT_SOURCE_DIR}/include")

# openmp support
find_package(OpenMP REQUIRED COMPONENTS C)
target_link_libraries(rbr OpenMP::OpenMP_C)
target_link_libraries(RBR OpenMP::OpenMP_C)

# blas support
find_package(BLAS REQUIRED)

# check whether it is MKL
if (BLAS_LIBRARIES MATCHES "libmkl_core")
    message(STATUS "Enable MKL support in RBR.")
    add_definitions(-DHAVE_MKL)
endif()

# RBR use cblas interface, so check whether it is included in BLAS
include(CheckFunctionExists)
set(CMAKE_REQUIRED_LIBRARIES "${BLAS_LIBRARIES}")
check_function_exists(cblas_daxpy BLAS_HAS_CBLAS)
target_link_libraries(RBR ${BLAS_LIBRARIES})
target_link_libraries(rbr ${BLAS_LIBRARIES})
# if not included, then try to find external cblas
if (BLAS_HAS_CBLAS)
    message(STATUS "BLAS library includes CBLAS.")
else()
    message(STATUS "BLAS library does not include CBLAS. Trying standalone CBLAS.")
    find_package(CBLAS REQUIRED)
    target_link_libraries(RBR CBLAS::CBLAS)
    target_link_libraries(rbr CBLAS::CBLAS)
endif()

# matlab support
if (BUILD_MATLAB_INTERFACE)
    find_package(Matlab COMPONENTS MEX_COMPILER)

    if (Matlab_FOUND)
        message(STATUS "Matlab found at ${Matlab_ROOT_DIR}")

        # add matlab interface library
        matlab_add_mex(
            NAME mex_rbr SHARED
            SRC ${MATLAB_SOURCE}
            LINK_TO RBR
        )

    else()
        message(WARNING "-DBUILD_MATLAB_INTERFACE=ON is set but Matlab not found. "
            "Try passing -DMatlab_ROOT_DIR to specify Matlab installation location")
    endif()
endif()


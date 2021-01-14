# - Find the standalone CBLAS library
#
#   Copyright (c) 2020, Haoyang Liu
#
# Usage:
#   find_package(CBLAS [REQUIRED] [QUIET])
#
# It sets the following variables:
#   CBLAS_FOUND                  ... true if CBLAS is found on the system
#   CBLAS_LIBRARIES              ... full paths to all found cblas libraries
#   CBLAS_INCLUDE_DIRS           ... full paths to cblas header files
#
# The following variables will be checked by the function
#   CBLAS_ROOT_DIR               ... if set, the libraries are exclusively searched
#                                    under this path
#

# possible install locations
set(CBLAS_HINTS
    ${CBLAS_ROOT_DIR}
    $ENV{CBLAS_DIR}
)
set(CBLAS_PATHS
    /usr
    /usr/local
    /System/Library/Framworks
)

# search header
find_path(CBLAS_INCLUDE_DIRS
    NAMES cblas.h
    HINTS ${CBLAS_HINTS}
    PATH_SUFFIXES
    include inc include/x86_64 include/x64
    openblas/include
    PATHS ${CBLAS_PATHS}
)
mark_as_advanced(CBLAS_INCLUDE_DIRS)

# search library
find_library(
    CBLAS_LIBRARIES
    NAMES cblas
    HINTS ${CBLAS_HINTS}
    PATH_SUFFIXES
    lib lib64 lib/x86_64 lib/x86
    openblas/lib
    PATHS ${CBlAS_PATHS}
)
mark_as_advanced(CBLAS_LIBRARIES)

if (NOT CBLAS_INCLUDE_DIRS)
    message(STATUS "Could NOT find CBLAS header file.")
endif()
if (NOT CBLAS_LIBRARIES)
    message(STATUS "Could NOT find CBLAS libraries.")
endif()

# define CBLAS interface target
if (CBLAS_INCLUDE_DIRS AND CBLAS_LIBRARIES)
    add_library(CBLAS::CBLAS INTERFACE IMPORTED)
    set_target_properties(CBLAS::CBLAS
        PROPERTIES INTERFACE_LINK_LIBRARIES "${CBLAS_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${CBLAS_INCLUDE_DIRS}"
        )
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CBLAS
    REQUIRED_VARS CBLAS_LIBRARIES CBLAS_INCLUDE_DIRS
    )


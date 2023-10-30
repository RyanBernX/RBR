# matlab support
if (BUILD_MATLAB_INTERFACE)
    find_package(Matlab COMPONENTS MEX_COMPILER)

    if (Matlab_FOUND)
        message(STATUS "Matlab found at ${Matlab_ROOT_DIR}")

        # add matlab interface library
        matlab_add_mex(
            NAME mex_rbr SHARED
            SRC mex_rbr.c
            LINK_TO RBR
        )

        # installing
        file(GLOB M_SCRIPTS *.m)
        install(TARGETS mex_rbr LIBRARY DESTINATION matlab RUNTIME DESTINATION matlab)
        install(FILES ${M_SCRIPTS} DESTINATION matlab)

    else()
        message(WARNING "-DBUILD_MATLAB_INTERFACE=ON is set but Matlab not found. "
            "Try passing -DMatlab_ROOT_DIR to specify Matlab installation location")
    endif()
endif()


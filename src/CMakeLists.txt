file(GLOB C_SOURCE *.c)
file(GLOB CXX_SOURCE *.cpp)

# main library
add_library(RBR ${C_SOURCE} ${CXX_SOURCE})
target_link_libraries(RBR PRIVATE ${CBLAS_LIBRARIES} OpenMP::OpenMP_C)

# RBR export headers
target_include_directories(RBR PUBLIC 
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
    "$<INSTALL_INTERFACE:include/RBR>")

# export symbols
if (WIN32)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wl,--exclude-all-symbols")
else()
    set_target_properties(RBR PROPERTIES
        C_VISIBILITY_PRESET hidden
        CXX_VISIBILITY_PRESET hidden
        VISIBILITY_INLINES_HIDDEN ON)
endif()


# subdirs
add_subdirectory(main)
add_subdirectory(matlab)

# installing
install(TARGETS RBR
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
install(FILES RBR.h adj_matrix.h api_macro.h types.h DESTINATION include/RBR)

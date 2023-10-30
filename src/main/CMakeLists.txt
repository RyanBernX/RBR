# define targets
add_executable(rbr main.c)
target_link_libraries(rbr PRIVATE RBR)

# runtime library path
if (APPLE)
    set_target_properties(rbr PROPERTIES INSTALL_RPATH @loader_path/../lib)
elseif (UNIX)
    set_target_properties(rbr PROPERTIES INSTALL_RPATH $ORIGIN/../lib)
endif()

# install 
install(TARGETS rbr RUNTIME DESTINATION bin)


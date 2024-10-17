#LIMA


# Now compile rest
add_executable(lima
    lima.cpp
 "render.h" "makebox.h" "argparser.h" "insertmolecules.h" "editconf.h")

target_include_directories(lima PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(lima LIMA_MD)


#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
#install(TARGETS lima DESTINATION ${CMAKE_SOURCE_DIR}/applications)

# Specify the installation target
if (UNIX)
    install(TARGETS lima DESTINATION /usr/bin)
else()
    install( TARGETS lima DESTINATION ${CMAKE_SOURCE_DIR}/applications )
endif()
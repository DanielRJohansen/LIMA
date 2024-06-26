#LIMA


# Now compile rest
add_executable(lima
    lima.cpp
 "selftest.h")

target_include_directories(lima PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(lima LIMA_MD)


#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
#install(TARGETS lima DESTINATION ${CMAKE_SOURCE_DIR}/applications)

# Specify the installation target
install(TARGETS lima DESTINATION /usr/bin)

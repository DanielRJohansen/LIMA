#LIMA_SERVER


# Now compile rest
add_executable(limaserver
    lima.cpp
)

target_include_directories(limaserver PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(limaserver LIMA_MD)


#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS limaserver DESTINATION ${CMAKE_SOURCE_DIR}/applications)
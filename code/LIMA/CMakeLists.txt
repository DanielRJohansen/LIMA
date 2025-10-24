#LIMA


# Now compile rest
add_executable(lima
    lima.cpp
  "include/SetWindowsRegistry.h")

target_include_directories(lima PRIVATE
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(lima LIMA_MD)


#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
#install(TARGETS lima DESTINATION ${CMAKE_SOURCE_DIR}/applications)

# Specify the installation target
if (UNIX)
    install(TARGETS lima DESTINATION /usr/bin)
else()
    install(TARGETS lima DESTINATION ${CMAKE_BINARY_DIR})

    target_sources(lima PRIVATE ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc)
endif()
#LIMA


# Now compile rest
add_executable(lima
    lima.cpp
    Daemon.cpp
  "include/WindowsIntegration.h")

target_include_directories(lima PRIVATE
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(lima LIMA_MD)
if (WIN32)
    target_link_libraries(lima ws2_32)
    # A running lima daemon locks lima.exe, so stop it before the linker overwrites the file
    add_custom_command(TARGET lima PRE_LINK
        COMMAND ${CMAKE_COMMAND} -DLIMA_EXE=$<TARGET_FILE:lima> -P ${CMAKE_CURRENT_SOURCE_DIR}/StopDaemon.cmake
        VERBATIM)
endif()
target_compile_definitions(lima PRIVATE LIMA_VERSION="${PROJECT_VERSION}")


#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
#install(TARGETS lima DESTINATION ${CMAKE_SOURCE_DIR}/applications)

# Specify the installation target
if (UNIX)
    install(TARGETS lima DESTINATION /usr/bin)
else()
    install(TARGETS lima DESTINATION ${CMAKE_BINARY_DIR})

    target_sources(lima PRIVATE ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc)
endif()

lima_pch(lima)

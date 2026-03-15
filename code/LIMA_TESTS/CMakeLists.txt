#limatest

if(WIN32)
    enable_language(RC)
    set(RESOURCE_FILES
        ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc
    )
endif()

set(SOURCE_FILES
    #src/AlgorithmTests.cu
    src/main.cpp
)
add_executable(limatest
    ${SOURCE_FILES}
    ${RESOURCE_FILES}
)

target_include_directories(limatest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(limatest LIMA_BASE)
target_link_libraries(limatest LIMA_MD)


if (UNIX)
    target_compile_options(limatest PRIVATE -O3)
    install(TARGETS limatest DESTINATION /usr/bin)
elseif (WIN32)
    #target_compile_options(LIMA_TESTS PRIVATE /O2)
endif()

message(STATUS "\n Resource files: ${RESOURCE_FILES}")


# Add the .exe for NCU profiling
add_executable(limaprofile
    src/limaprofile.cpp
)
target_include_directories(limaprofile PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(limaprofile LIMA_BASE LIMA_MD)
#if (UNIX OR WIN32)
#    target_compile_options(limaprofile PRIVATE -O3)
#    install(TARGETS limaprofile DESTINATION /usr/bin)
#elseif (WIN32)
#    #target_compile_options(limaprofile PRIVATE /O2)
#endif()

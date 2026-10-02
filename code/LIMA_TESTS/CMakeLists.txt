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

option(ENABLE_EXPERIMENTAL_BATCH_TESTS "Enable tests that require shared-kernel batch execution" OFF)
if (ENABLE_EXPERIMENTAL_BATCH_TESTS)
    target_compile_definitions(limatest PRIVATE ENABLE_EXPERIMENTAL_BATCH_TESTS)
endif()

target_include_directories(limatest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(limatest LIMA_BASE)
target_link_libraries(limatest LIMA_MD)

add_executable(agenttesting
    src/agenttesting.cpp
    src/EnergyMinimizationtests.cpp
    ${RESOURCE_FILES}
)
target_include_directories(agenttesting PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(agenttesting LIMA_BASE LIMA_MD)

# Release tests for the lima CLI programs. Runs lima.exe as a child process, so Windows only for now
if (WIN32)
    add_executable(limaclitest
        src/CliTests.cpp
        ${RESOURCE_FILES}
    )
    target_include_directories(limaclitest PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/code/LIMA/include
    )
    target_link_libraries(limaclitest LIMA_BASE LIMA_MD)
    add_dependencies(limaclitest lima)
endif()


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

lima_pch(limatest)
lima_pch(agenttesting)
lima_pch(limaprofile)

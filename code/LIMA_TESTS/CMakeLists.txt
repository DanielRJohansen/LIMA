#limatest

if(WIN32)
    enable_language(RC)
    set(RESOURCE_FILES
        ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc
    )
endif()

# The tests, compiled once and shared by the test executables. Declared in Tests.h
add_library(LIMA_TESTCASES STATIC
    src/AlgorithmTests.cpp
    src/BatchingTests.cpp
    src/Benchmarks.cpp
    src/ElectrostaticsTests.cpp
    src/EngineBatchTests.cpp
    src/FileTests.cpp
    src/ForceComparisons.cpp
    src/ForceCorrectness.cpp
    src/ForcefieldTests.cpp
    src/MDStability.cpp
    src/MembraneBuilder.cpp
    src/MinorPrograms.cpp
    src/ProgramsTests.cpp
    src/SetupTests.cpp
)
target_include_directories(LIMA_TESTCASES PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(LIMA_TESTCASES PUBLIC LIMA_BASE LIMA_MD)
# The test sources are compiled as a few unity batches (see lima_unity), so TestUtils.h is parsed a few times only, and
# a test-specific PCH would cost more to build than it saves
lima_pch(LIMA_TESTCASES)

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

target_link_libraries(limatest LIMA_TESTCASES)

add_executable(agenttesting
    src/agenttesting.cpp
    src/EnergyMinimizationtests.cpp
    ${RESOURCE_FILES}
)
target_include_directories(agenttesting PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_link_libraries(agenttesting LIMA_TESTCASES)

# Release tests for the lima CLI programs. Runs lima.exe as a child process, so Windows only for now
if (WIN32)
    # Compiled as an object library, so only linking waits for lima.exe (add_dependencies below), not compiling
    add_library(limaclitest_objects OBJECT src/CliTests.cpp)
    target_include_directories(limaclitest_objects PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_SOURCE_DIR}/code/LIMA/include
    )
    target_link_libraries(limaclitest_objects PUBLIC LIMA_BASE LIMA_MD)
    add_executable(limaclitest
        $<TARGET_OBJECTS:limaclitest_objects>
        ${RESOURCE_FILES}
    )
    target_link_libraries(limaclitest limaclitest_objects)
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
target_link_libraries(limaprofile LIMA_TESTCASES)
#if (UNIX OR WIN32)
#    target_compile_options(limaprofile PRIVATE -O3)
#    install(TARGETS limaprofile DESTINATION /usr/bin)
#elseif (WIN32)
#    #target_compile_options(limaprofile PRIVATE /O2)
#endif()

lima_pch(limatest)
lima_pch(agenttesting)
lima_pch(limaprofile)

#LIMA_TESTS
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

#file(GLOB LIMA_TESTS_SOURCES "src/*.cpp")

add_executable(limatests
    src/main.cpp
    src/EnergyMinimization.cuh
    src/ForceCorrectness.h
    src/MDStability.cuh
    src/TestUtils.h
)

target_include_directories(limatests PRIVATE $ {CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(limatests LIMA_BASE)

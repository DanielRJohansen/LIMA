#LIMA_ENGINE

find_package(CUDA REQUIRED)

add_library(LIMA_ENGINE
    src/Engine.cu
    src/EngineKernels.cu
    src/EngineUtils.cu
    src/Neighborlists.cu
    src/Thermostat.cpp
)

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
${CUDA_INCLUDE_DIRS}
)

target_link_libraries(LIMA_ENGINE LIMA_BASE)

set_target_properties(LIMA_ENGINE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_ENGINE PROPERTIES CUDA_ARCHITECTURES "89")

# Self-test executable
add_executable(engine_self_test src/Main.cpp)
target_link_libraries(engine_self_test LIMA_ENGINE)

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
${CUDA_INCLUDE_DIRS}
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_ENGINE LIMA_BASE)

set_target_properties(LIMA_ENGINE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_ENGINE PROPERTIES CUDA_ARCHITECTURES "89")

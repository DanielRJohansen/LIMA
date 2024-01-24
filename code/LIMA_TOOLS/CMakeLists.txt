#LIMA_TOOLS

find_package(CUDA REQUIRED)

set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_89")


add_library(LIMA_TOOLS
    src/SimulationBuilder.cpp
    src/VirtualPathMaker.cpp
    src/CompoundBuilder.cpp
    )

target_include_directories(LIMA_TOOLS PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_TOOLS PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
${CUDA_INCLUDE_DIRS}
)

target_link_libraries(LIMA_TOOLS LIMA_BASE)
target_link_libraries(LIMA_TOOLS LIMA_FORCEFIELDMAKER)

set_target_properties(LIMA_TOOLS PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_TOOLS PROPERTIES CUDA_ARCHITECTURES "89")

if (UNIX)
    target_compile_options(LIMA_TOOLS PRIVATE -use_fast_math -O3)
endif()



# Self-test executable
# add_executable(tools_self_test src/EngineMain.cpp "src/EngineBodies.cuh" "src/BoundaryConditionPublic.cpp")
# target_link_libraries(tools_self_test LIMA_TOOLS)

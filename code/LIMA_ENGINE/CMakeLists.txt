#LIMA_ENGINE
add_library(LIMA_ENGINE
    src/Engine.cu
    src/EngineKernels.cu
    src/EngineUtils.cu
    src/Neighborlists.cu 
    src/Thermostat.cu
    src/BoundaryConditionPublic.cpp)

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
${CUDA_INCLUDE_DIRS}
)

target_link_libraries(LIMA_ENGINE LIMA_BASE)

if (UNIX)
    target_compile_options(LIMA_ENGINE PRIVATE -use_fast_math -O3)
endif()



# Self-test executable
add_executable(engine_self_test src/EngineMain.cu)
target_link_libraries(engine_self_test LIMA_ENGINE)

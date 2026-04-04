#LIMA_ENGINE

# Do NOT add any other .cu files, the engine should be a single compilation unit
# in order to access __constant__ memory
add_library(LIMA_ENGINE STATIC
    src/Engine.cu
    src/SimulationDevice.cpp
  "src/ParticleClusteringCPU.cpp"
)

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_ENGINE LIMA_BASE CUDA::cufft)


if (UNIX)
    target_compile_options(LIMA_ENGINE PRIVATE -use_fast_math -O3)
elseif (WIN32)
       target_compile_options(LIMA_ENGINE PRIVATE 
        $<$<CONFIG:Release>:/use_fast_math>
        $<$<CONFIG:Release>:/O2>
        $<$<COMPILE_LANGUAGE:CUDA>:--maxrregcount=255>
    )
endif()

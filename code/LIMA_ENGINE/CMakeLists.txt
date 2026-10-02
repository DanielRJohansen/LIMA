#LIMA_ENGINE

# Do NOT add any other .cu files with kernels, the engine should be a single compilation unit in order to access
# __constant__ memory. ThrustWrappers.cu is the exception: it has no __constant__ access, and keeping thrust's
# headers out of Engine.cu makes that (the slowest file in the build) compile ~16 s faster
add_library(LIMA_ENGINE STATIC
    src/Engine.cu
    src/ThrustWrappers.cu
    src/SimulationData.cpp
    src/EnergyMinimizationHost.cpp
  "src/ParticleClusteringCPU.cpp"
)

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

# Static on Linux, so the packages only need the NVIDIA driver and no CUDA libraries.
# Windows has no static cuFFT; its DLL is shipped next to lima.exe
if (UNIX)
    target_link_libraries(LIMA_ENGINE LIMA_BASE CUDA::cufft_static_nocallback)	# The callback variant needs device linking, and LIMA uses no cuFFT callbacks
else()
    target_link_libraries(LIMA_ENGINE LIMA_BASE CUDA::cufft)
endif()


if (UNIX)
    target_compile_options(LIMA_ENGINE PRIVATE -use_fast_math -O3)
elseif (WIN32)
       target_compile_options(LIMA_ENGINE PRIVATE 
        $<$<COMPILE_LANGUAGE:CUDA>:-use_fast_math>
        $<$<AND:$<CONFIG:Release>,$<COMPILE_LANGUAGE:CXX>>:/O2>
        $<$<COMPILE_LANGUAGE:CUDA>:--maxrregcount=255>
    )
endif()

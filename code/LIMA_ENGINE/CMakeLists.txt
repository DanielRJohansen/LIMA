#LIMA_ENGINE

# The engine is split in several .cu files so they compile in parallel, as Engine.cu is the longest step of the build.
# This requires that no __constant__ or other __device__ variables are shared between them, as there is no
# relocatable device code. Each kernel must be defined in one .cu file only.
add_library(LIMA_ENGINE STATIC
    src/Engine.cu
    src/PME.cu
    src/ParticleClusters.cu
    src/SuperclusterTaskBuilder.cu
    src/CubWrappers.cu
    src/SimulationData.cpp
    src/EngineBatch.cpp
    src/EnergyMinimizationHost.cpp
  "src/ParticleClusteringCPU.cpp"
)

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (LIMA_GPU_BACKEND STREQUAL "HIP")
    # The same .cu sources, compiled as HIP
    get_target_property(engineSources LIMA_ENGINE SOURCES)
    list(FILTER engineSources INCLUDE REGEX "\\.cu$")
    set_source_files_properties(${engineSources} PROPERTIES LANGUAGE HIP)
    find_package(hipfft REQUIRED CONFIG)
    find_package(hipcub REQUIRED CONFIG)
    target_link_libraries(LIMA_ENGINE LIMA_BASE hip::hipfft hip::hipcub)
    # The equivalent of the CUDA build's -use_fast_math: flush denormals, approximate fp32 division/sqrt and math functions.
    # Not -ffast-math, whose finite-math assumption would remove the kernels' NaN checks. T4 on the 780M: +16%
    set(hipFastMath -fgpu-flush-denormals-to-zero -fno-hip-fp32-correctly-rounded-divide-sqrt -fapprox-func)
    target_compile_options(LIMA_ENGINE PRIVATE "$<$<COMPILE_LANGUAGE:HIP>:${hipFastMath}>")
    return()	# The rest is the CUDA build's flags
endif()

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

# Main cmakelist
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(SET CMP0146 NEW) # Set the policy to NEW

if(UNIX)
    set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
endif()

# Dev builds target only the dev GPU (RTX 40-series) to keep compile times down. Release builds pass
# -DCMAKE_CUDA_ARCHITECTURES=... to target every supported GPU, see distribution/release.bat
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES "89")
endif()

# CUDA (NVIDIA) or HIP (AMD). Both compile the same sources; GpuRuntime.h maps the CUDA runtime names to HIP. HIP needs the
# ROCm SDK (-DCMAKE_HIP_COMPILER=<ROCm>/bin/clang++.exe -DCMAKE_PREFIX_PATH=<ROCm>), see the x64-Release-HIP configuration
set(LIMA_GPU_BACKEND "CUDA" CACHE STRING "GPU backend: CUDA or HIP")
set_property(CACHE LIMA_GPU_BACKEND PROPERTY STRINGS CUDA HIP)
if (NOT LIMA_GPU_BACKEND MATCHES "^(CUDA|HIP)$")
	message(FATAL_ERROR "LIMA_GPU_BACKEND must be CUDA or HIP, was ${LIMA_GPU_BACKEND}")
endif()
if (LIMA_GPU_BACKEND STREQUAL "HIP" AND NOT DEFINED CMAKE_HIP_ARCHITECTURES)
	set(CMAKE_HIP_ARCHITECTURES "gfx1103")	# The dev AMD GPU (Radeon 780M)
endif()
set(GPU_TARGETS ${CMAKE_HIP_ARCHITECTURES})	# For ROCm's package configs, which otherwise probe for a GPU

# The single source of the release version. distribution/release.bat publishes whatever is set here
project(lima VERSION 1.4.0 LANGUAGES CXX ${LIMA_GPU_BACKEND} C)

# CMake can't mix Clang and MSVC, so HIP builds compile all C/C++ with ROCm's clang (on the MSVC STL and ABI)
if (LIMA_GPU_BACKEND STREQUAL "HIP" AND CMAKE_CXX_SIMULATE_ID STREQUAL "MSVC" AND CMAKE_CXX_SIMULATE_VERSION VERSION_GREATER_EQUAL 19.50)
	message(FATAL_ERROR "ROCm 7.2's clang can't parse the MSVC 14.5x STL: configure from a VS 2022 (14.44) developer shell")
endif()

if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:preprocessor /Zc:__cplusplus")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Zc:preprocessor")

	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=/Zc:preprocessor -Xcompiler=/Zc:__cplusplus")
endif()

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 23)
# Not 23: from C++23 the MSVC STL makes isinf/isnan/... constexpr, which HIP treats as host+device, clashing with its own
set(CMAKE_HIP_STANDARD 20)
#set(CMAKE_BUILD_TYPE Debug)

#set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CUDA_SEPARABLE_COMPILATION ON)


set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_89,code=sm_89")
# 20012: glm defaulted-ctor __host__/__device__ annotations (third-party); 20054: dynamic init of __shared__ variables (intended)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr -lineinfo -diag-suppress=20012,20054") # --ptxas-options=-v")


if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++26 -static-libgcc -static-libstdc++" )
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,relro,-z,now")	#Flags for executable security, to appease NAMCAP
else()
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std=c++26 /W4" )
	#add_compile_options( $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/Zc:preprocessor>)
	#add_compile_options(/std:c++latest /Zc:__cplusplus /W4)
    #add_compile_options(
    #    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=/Zc:preprocessor>
    #    $<$<COMPILE_LANG_AND_ID:CXX,MSVC>:/Zc:preprocessor>
    #)
endif()


find_package(OpenGL REQUIRED)
if (LIMA_GPU_BACKEND STREQUAL "HIP")
	find_package(hip REQUIRED CONFIG)
	# hip::device (linked through hipFFT/hipCUB) would compile every C++ file as HIP (-x hip). CMake's HIP language already
	# compiles the .cu files as HIP for CMAKE_HIP_ARCHITECTURES, so the plain C++ files stay C++
	set_property(TARGET hip::device PROPERTY INTERFACE_COMPILE_OPTIONS "")
	add_compile_definitions(LIMA_HIP)
	set(LIMA_GPU_INCLUDE_DIRS "")			# hip::host, linked by LIMA_BASE, brings the HIP headers
else()
	find_package(CUDAToolkit REQUIRED)
	set(LIMA_GPU_INCLUDE_DIRS ${CUDAToolkit_INCLUDE_DIRS})
endif()
#find_package(PythonLibs REQUIRED)

if (UNIX)
    find_package(TBB REQUIRED)
    find_package(glm REQUIRED)
    find_package(glfw3 REQUIRED)
endif()

# Set the installation prefix
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/applications")
 
add_compile_definitions(GLEW_STATIC)
# Thrust/CUB otherwise emit NVTX profiler ranges, which drags windows.h into every file using them (~23% of Engine.cu's parsing)
add_compile_definitions(CCCL_DISABLE_NVTX)

# Lets a target's objects compile without waiting for the static libraries it links to
set(CMAKE_OPTIMIZE_DEPENDENCIES ON)

# Precompiled header for the host (C++) code of a target. CUDA files are left alone, nvcc handles PCH poorly.
# MSVC lowercases the .pch path, which fails if the build tree is inside a case-sensitive directory
# (fsutil file queryCaseSensitiveInfo <dir>). Configure with -DLIMA_PCH=OFF there
option(LIMA_PCH "Use a precompiled header for host code" ON)
function(lima_pch target)
	if (LIMA_PCH)
		if (target STREQUAL "LIMA_BASE")
			target_precompile_headers(LIMA_BASE PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${CMAKE_SOURCE_DIR}/code/LIMA_BASE/include/lima_pch.h>")
		else()
			# One PCH for the whole project, instead of ~7 s to build one per target
			target_precompile_headers(${target} REUSE_FROM LIMA_BASE)
		endif()
	endif()
endfunction()

# Unity build: a target's C++ sources are compiled in batches, which saves parsing the same headers and compiling the
# same template instantiations once per file. The CUDA targets are excluded, their sources must stay separate.
# A batch should not take much longer than Engine.cu, the longest step of the build
option(LIMA_UNITY "Compile the C++ sources of a target in batches" ON)
function(lima_unity target batchSize)
	if (LIMA_UNITY)
		set_target_properties(${target} PROPERTIES UNITY_BUILD ON UNITY_BUILD_BATCH_SIZE ${batchSize})
	endif()
endfunction()

# Ninja starts jobs in the order targets are defined, so the CUDA libraries go first: Engine.cu is the slowest
# compile, and every other target can be built alongside it
add_subdirectory(code/LIMA_ENGINE)
add_subdirectory(code/LIMA_CONVEXHULLENGINE)
add_subdirectory(code/LIMA_DEPENDS)
add_subdirectory(code/LIMA_BASE)
add_subdirectory(code/LIMA_TOOLS)
add_subdirectory(code/LIMA_MD)
add_subdirectory(code/LIMA_TESTS)
add_subdirectory(code/LIMA)

# These include <windows.h> on Windows, whose macros would leak into the sources batched after them
if (WIN32)
	set_source_files_properties(code/LIMA_BASE/src/Filehandling.cpp DIRECTORY code/LIMA_BASE PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
	set_source_files_properties(code/LIMA_MD/src/Analyzer.cpp DIRECTORY code/LIMA_MD PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
endif()
# The longest C++ sources compile on their own, so no batch gets much longer than them
set_source_files_properties(code/LIMA_MD/src/Pdb2gmx.cpp code/LIMA_MD/display/src/Display.cpp DIRECTORY code/LIMA_MD PROPERTIES SKIP_UNITY_BUILD_INCLUSION ON)
lima_unity(LIMA_BASE 8)
lima_unity(LIMA_TOOLS 8)
lima_unity(LIMA_MD 5)
lima_unity(LIMA_TESTCASES 5)
lima_unity(agenttesting 8)

if (UNIX)
	# Copy the resources to the appropriate directory
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/ DESTINATION /usr/share/LIMA/resources)
	#install(DIRECTORY ${CMAKE_SOURCE_DIR}/dev/ DESTINATION /usr/share/LIMA/dev)
else()
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/ DESTINATION ${CMAKE_BINARY_DIR}/resources)
	#install(DIRECTORY ${CMAKE_SOURCE_DIR}/dev/ DESTINATION ${CMAKE_BINARY_DIR}/dev)

	# Run the program once after install to register file association
    install(CODE "
        execute_process(
            COMMAND \"\$<TARGET_FILE:lima>\" setregistry
            RESULT_VARIABLE result
        )
        if(NOT result EQUAL 0)
            message(WARNING \"LIMA registry setup failed with code \${result}\")
        endif()
    ")
endif()

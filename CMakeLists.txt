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

# The single source of the release version. distribution/release.bat publishes whatever is set here
project(lima VERSION 1.3.0 LANGUAGES CXX CUDA C)

if (MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /Zc:preprocessor /Zc:__cplusplus")
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} /Zc:preprocessor")

	set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=/Zc:preprocessor -Xcompiler=/Zc:__cplusplus")
endif()

set(CMAKE_CXX_STANDARD 26)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 23)
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
find_package(CUDAToolkit REQUIRED)
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

# Ninja starts jobs in the order targets are defined, so the CUDA libraries go first: they are the slowest
# compiles by far (~70 s each), and every other target can be built alongside them
add_subdirectory(code/LIMA_ENGINE)
add_subdirectory(code/LIMA_CONVEXHULLENGINE)
add_subdirectory(code/LIMA_DEPENDS)
add_subdirectory(code/LIMA_BASE)
add_subdirectory(code/LIMA_TOOLS)
add_subdirectory(code/LIMA_MD)
add_subdirectory(code/LIMA_TESTS)
add_subdirectory(code/LIMA)

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

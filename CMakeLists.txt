# Main cmakelist
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(SET CMP0146 NEW) # Set the policy to NEW

if(UNIX)
    set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
endif()

set(CMAKE_CUDA_ARCHITECTURES "89" "75" "61")

project(lima LANGUAGES CXX CUDA C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)


set(CUDA_SEPARABLE_COMPILATION ON)


set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_89,code=sm_89")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr")

# Only gcc v13 and newer can compile LIMA. On arch gcc defaults to v13, but on ubuntu we need to explicitly install it.
#find_program(GCC_13_FOUND g++-13)
#if(GCC_13_FOUND)
#	set(CMAKE_CXX_COMPILER g++-13)
#endif()

if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -static-libgcc -static-libstdc++" )
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,relro,-z,now")	#Flags for executable security, to appease NAMCAP
else()
	set(set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std=c++20 /W4" ))
endif()


find_package(OpenGL REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(PythonLibs REQUIRED)

#include_directories(${PYTHON_INCLUDE_DIRS})
#link_directories(${PYTHON_LIBRARIES})

if (UNIX)
	find_package(glm REQUIRED)
	find_package(glfw3 REQUIRED)
endif()

# Set the installation prefix
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/applications")

add_compile_definitions(GLEW_STATIC)

add_subdirectory(code/LIMA_DEPENDS)
add_subdirectory(code/LIMA_BASE)
add_subdirectory(code/LIMA_ENGINE)
add_subdirectory(code/LIMA_TOOLS)
add_subdirectory(code/LIMA_MD)
add_subdirectory(code/LIMA_TESTS)
add_subdirectory(code/LIMA)

if (UNIX)
	# Copy the resources to the appropriate directory
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/ DESTINATION /opt/LIMA/resources)

	# Add a post-installation step to set the permissions
	install(CODE "execute_process(COMMAND chmod -R a+rw /opt/LIMA/resources)")
endif()

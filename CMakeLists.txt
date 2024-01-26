# Main cmakelist
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

if(UNIX)
    set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
endif()

set(CMAKE_CUDA_ARCHITECTURES "89" "75" "61")

project(lima LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)


set(CUDA_SEPARABLE_COMPILATION ON)


#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -arch=sm_89")
set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_61,code=sm_61 -gencode arch=compute_75,code=sm_75 -gencode arch=compute_89,code=sm_89")

if (UNIX)
	set(CMAKE_CXX_COMPILER g++-13)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
endif()

find_package(CUDA REQUIRED)
find_package(OpenGL REQUIRED)

# Set the installation prefix
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../applications")

add_subdirectory(code/LIMA_BASE)
add_subdirectory(code/LIMA_ENGINE)
add_subdirectory(code/LIMA_TOOLS)
add_subdirectory(code/LIMA_FORCEFIELDMAKER)
add_subdirectory(code/LIMA_MD)
add_subdirectory(code/LIMA_TESTS)
add_subdirectory(code/LIMA_SERVER)



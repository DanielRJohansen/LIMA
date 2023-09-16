# Main cmakelist
cmake_minimum_required(VERSION 3.26 FATAL_ERROR)

if(UNIX)
    set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
endif()

set(CMAKE_CUDA_ARCHITECTURES "89")

project(lima LANGUAGES CXX CUDA)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)

if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
endif()
find_package(OpenGL REQUIRED)

# Set the installation prefix
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/../applications")

add_subdirectory(LIMA_BASE)
add_subdirectory(LIMA_ENGINE)
add_subdirectory(LIMA_FORCEFIELDMAKER)
add_subdirectory(LIMA_MD)
add_subdirectory(LIMA_TESTS)
add_subdirectory(LIMA_SERVER)



#LIMA_APPS
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

# We dont need CUDA, but we need the CMAKE error as soon as possible
project(lima LANGUAGES CXX CUDA)

# Set the installation directory to /usr/local/bin
set(CMAKE_INSTALL_PREFIX /usr/local)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)

if (UNIX)
    set(CMAKE_CXX_COMPILER g++-13)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
endif()


add_executable(lima
    lima.cpp
)

target_include_directories(lima PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# Specify the installation target
install(TARGETS lima DESTINATION bin)

#LIMA_APPS
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

project(lima LANGUAGES CXX)

# Set the installation directory to /usr/local/bin
set(CMAKE_INSTALL_PREFIX /usr/local)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)

if (UNIX)
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
install(TARGETS lima DESTINATION /bin)

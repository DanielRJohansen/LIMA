# Main cmakelist
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)
cmake_policy(SET CMP0146 NEW) # Set the policy to NEW

if(UNIX)
    set(CMAKE_INSTALL_MESSAGE NEVER)  # Install silently.
endif()

set(CMAKE_CUDA_ARCHITECTURES "89")

project(lima LANGUAGES CXX CUDA C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
#set(CMAKE_BUILD_TYPE Debug)


set(CUDA_SEPARABLE_COMPILATION ON)


set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -gencode arch=compute_89,code=sm_89")
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --expt-relaxed-constexpr") # --ptxas-options=-v")


if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20 -static-libgcc -static-libstdc++" )
	add_definitions(-w)	# Hides the shared varialbe dynamic initilization error in kernels
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,-z,relro,-z,now")	#Flags for executable security, to appease NAMCAP
else()
	set(set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std=c++20 /W4" ))
	#set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "${CMAKE_COMMAND} -E time")
	#add_compile_options(/external:W0 /external:anglebrackets)
endif()



find_package(OpenGL REQUIRED)
find_package(CUDAToolkit REQUIRED)
find_package(PythonLibs REQUIRED)

if (UNIX)
	find_package(TBB REQUIRED)
	find_package(glm REQUIRED)
	find_package(glfw3 REQUIRED)
endif()

# Set the installation prefix
set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_SOURCE_DIR}/applications")

add_compile_definitions(GLEW_STATIC)

add_subdirectory(code/LIMA_DEPENDS)
add_subdirectory(code/LIMA_BASE)
add_subdirectory(code/LIMA_ENGINE)
add_subdirectory(code/LIMA_CONVEXHULLENGINE)
add_subdirectory(code/LIMA_TOOLS)
add_subdirectory(code/LIMA_MD)
add_subdirectory(code/LIMA_TESTS)
add_subdirectory(code/LIMA)

if (UNIX)
	# Copy the resources to the appropriate directory
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/resources/ DESTINATION /usr/share/LIMA/resources)
	install(DIRECTORY ${CMAKE_SOURCE_DIR}/dev/ DESTINATION /usr/share/LIMA/dev)
endif()

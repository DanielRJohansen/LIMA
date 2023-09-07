# Main cmakelist
cmake_minimum_required(VERSION 3.26 FATAL_ERROR)

#set(CUDAToolkit_ROOT /opt/cuda/bin/)
set(CMAKE_CUDA_ARCHITECTURES "89")
#set(CMAKE_CUDA_COMPILER /opt/cuda/bin/nvcc)

project(lima LANGUAGES CXX CUDA)

#set_target_properties(lima PROPERTIES CUDA_ARCHITECTURES "89")


set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)
#set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -j4")  # Limit to 4 concurrent tasks

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20")



#find_package(CUDA REQUIRED)
#find_package(Chrono REQUIRED)


#add_subdirectory(${PROJECT_SOURCE_DIR}/dependencies/GLFW/)
find_package(OpenGL REQUIRED)


add_subdirectory(LIMA_BASE)
add_subdirectory(LIMA_ENGINE)
add_subdirectory(LIMA_FORCEFIELDMAKER)
add_subdirectory(LIMA_MD)
add_subdirectory(LIMA_TESTS)


#include_directories("${CUDA_INCLUDE_DIRS}")
#set(CUDA_LIBRARIES ${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcudart.so)
#target_link_libraries(lima glfw)
#target_link_libraries(lima OpenGL::GL)

#target_compile_features(lima PUBLIC cxx_std_20)
#set_target_properties(lima PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(lima PROPERTIES CUDA_ARCHITECTURES "89")


#target_include_directories(lima PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(lima PUBLIC "${CUDA_INCLUDE_DIRS}")

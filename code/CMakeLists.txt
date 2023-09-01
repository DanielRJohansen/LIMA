# Main cmakelist
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

set(CUDAToolkit_ROOT /opt/cuda/bin/)
#set(CUDA_ARCHITECTURES "89")
#set_target_properties(lima PROPERTIES CUDA_ARCHITECTURES "89")
set(CMAKE_CUDA_ARCHITECTURES 89)
set(CMAKE_CUDA_COMPILER /opt/cuda/bin/nvcc)

project(lima LANGUAGES CXX CUDA)



set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_REQUIRED ON)
set(CMAKE_CUDA_STANDARD 20)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++20")



#find_package(CUDA REQUIRED)
#find_package(Chrono)


#add_subdirectory(${PROJECT_SOURCE_DIR}/dependencies/GLFW/)
find_package(OpenGL REQUIRED)


add_subdirectory(src/LIMA_BASE)
#add_subdirectory(LIMA_ENGINE)
#add_subdirectory(LIMA_FORCFIELDMAKER)
#add_subdirectory(LIMA_MD)
add_subdirectory(src/LIMA_TESTS)


#target_link_libraries(lima glfw)
#target_link_libraries(lima OpenGL::GL)

#target_compile_features(lima PUBLIC cxx_std_20)
#set_target_properties(lima PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(lima PROPERTIES CUDA_ARCHITECTURES "61;75")


#target_include_directories(lima PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
#target_include_directories(lima PUBLIC "${CUDA_INCLUDE_DIRS}")

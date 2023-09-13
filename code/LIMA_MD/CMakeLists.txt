#LIMA_MD
find_package(CUDA REQUIRED)
find_package(OpenGL REQUIRED)

add_library(LIMA_MD STATIC
    src/Analyzer.cu
    src/BoxBuilder.cu
    src/ChemfilesInterface.cpp
    src/CompoundBuilder.cpp
    src/DisplayV2.cpp
    src/Environment.cpp
    src/Interface.cpp
    src/Rasterizer.cu
    src/VirtualPathMaker.cpp
)


target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}    
    ${CMAKE_SOURCE_DIR}/dependencies/GLFW/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_include_directories(LIMA_MD PRIVATE

)



target_link_libraries(LIMA_MD LIMA_BASE)
target_link_libraries(LIMA_MD LIMA_ENGINE)
target_link_libraries(LIMA_MD LIMA_FORCEFIELDMAKER)
if (WIN32)
    target_link_libraries(LIMA_MD ${CMAKE_SOURCE_DIR}/dependencies/GLFW/lib/glfw3.lib)
endif()
#target_include_directories(LIMA_MD PRIVATE ${OPENGL_INCLUDE_DIR})
target_link_libraries(LIMA_MD ${OPENGL_LIBRARIES})


set_target_properties(LIMA_MD PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_MD PROPERTIES CUDA_ARCHITECTURES "89")

if (UNIX)
    target_compile_options(LIMA_MD PRIVATE -use_fast_math -O3)
endif()





# Self-test executable
#set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --relocatable-device-code=true")
find_package(CUDA REQUIRED)
add_executable(md_self_test src/Main.cpp)
set_target_properties(md_self_test PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(md_self_test PROPERTIES CUDA_ARCHITECTURES "89")
#target_link_libraries(md_self_test LIMA_BASE)
#target_link_libraries(md_self_test LIMA_ENGINE)
#target_link_libraries(md_self_test LIMA_FORCEFIELDMAKER)
target_link_libraries(md_self_test LIMA_MD)

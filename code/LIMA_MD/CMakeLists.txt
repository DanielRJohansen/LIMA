#LIMA_MD
#find_package(OpenGL REQUIRED)

add_library(LIMA_MD STATIC
    src/Analyzer.cu
    src/BoxBuilder.cpp
    src/DisplayV2.cpp
    src/Environment.cpp
    src/Interface.cpp
    src/Rasterizer.cu
    
)


target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}    
    ${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)
target_include_directories(LIMA_MD PRIVATE
)



target_link_libraries(LIMA_MD LIMA_BASE)
target_link_libraries(LIMA_MD LIMA_ENGINE)
target_link_libraries(LIMA_MD LIMA_TOOLS)
target_link_libraries(LIMA_MD LIMA_FORCEFIELDMAKER)
if (WIN32)
    target_link_libraries(LIMA_MD ${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/lib/glfw3.lib)
endif()
#target_include_directories(LIMA_MD PRIVATE ${OPENGL_INCLUDE_DIR})
target_link_libraries(LIMA_MD ${OPENGL_LIBRARIES})



if (UNIX)
    target_compile_options(LIMA_MD PRIVATE -use_fast_math -O3)
endif()





# Self-test executable
add_executable(md_self_test src/Main.cpp)
target_link_libraries(md_self_test LIMA_MD)

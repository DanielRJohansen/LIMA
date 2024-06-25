#LIMA_MD

add_library(LIMA_MD STATIC
    src/Analyzer.cu
    src/BoxBuilder.cpp
    src/DisplayV2.cpp
    src/Environment.cpp
    src/Interface.cpp
    src/Rasterizer.cu
    src/Programs.cpp
    src/Render.cpp
)


# Set the path to GLEW installation
set(GLEW_LIBRARY_DIR "${CMAKE_SOURCE_DIR}/code/dependencies/glew-2.1.0/lib/Release/x64")
set(GLEW_LIBRARIES "${GLEW_LIBRARY_DIR}/glew32s.lib")


target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}    
    ${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/include
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/code/dependencies/glew-2.1.0/include
    ${CMAKE_SOURCE_DIR}/code/dependencies/glm/
)
target_include_directories(LIMA_MD PRIVATE
)



target_link_libraries(LIMA_MD LIMA_BASE)
target_link_libraries(LIMA_MD LIMA_ENGINE)
target_link_libraries(LIMA_MD LIMA_TOOLS)

target_link_libraries(LIMA_MD ${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/lib/glfw3.lib)

target_link_libraries(LIMA_MD ${OPENGL_LIBRARIES})

# Define GLEW_STATIC to indicate that you are using the static version of GLEW
target_compile_definitions(LIMA_MD PUBLIC GLEW_STATIC)

# Link GLEW libraries to the executable
target_link_libraries(LIMA_MD ${GLEW_LIBRARIES})


if (UNIX)
    target_compile_options(LIMA_MD PRIVATE -use_fast_math -O3)
endif()





# Self-test executable
add_executable(md_self_test src/Main.cpp)
target_link_libraries(md_self_test LIMA_MD)

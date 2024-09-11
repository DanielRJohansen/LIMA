#LIMA_MD

add_library(LIMA_MD STATIC
    src/Analyzer.cu
    src/BoxBuilder.cpp
    src/Display.cpp
    src/Environment.cpp
    src/Interface.cpp
    src/Programs.cpp
    src/Rendering.cpp
)



if (NOT UNIX)
    set(GLFW_LIBRARIES "${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/lib/glfw3.lib")
    set(GLFW_INCLUDE_DIRS "${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/include")
else()
    set(GLFW_LIBRARIES glfw)
endif()


target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${GLFW_INCLUDE_DIRS}
)
target_include_directories(LIMA_MD PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CMAKE_SOURCE_DIR}/code/dependencies/glm/
)


# Define GLEW_STATIC to indicate that you are using the static version of GLEW
target_compile_definitions(LIMA_MD PUBLIC GLEW_STATIC)


target_link_libraries(LIMA_MD 
    LIMA_BASE
    LIMA_ENGINE
    LIMA_TOOLS
    LIMA_DEPENDS
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
)


if (UNIX)
    target_compile_options(LIMA_MD PRIVATE -use_fast_math -O3)
endif()





# Self-test executable
add_executable(md_self_test src/Main.cpp "src/Shaders.h")
target_link_libraries(md_self_test LIMA_MD)

#LIMA_TESTS

find_package(CUDA REQUIRED)

add_executable(limatests
    src/main.cpp
)

target_include_directories(limatests PRIVATE
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)


target_link_libraries(limatests LIMA_BASE)
target_link_libraries(limatests LIMA_MD)

set_target_properties(limatests PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(limatests PROPERTIES CUDA_ARCHITECTURES "89")

if (UNIX)
    target_compile_options(limatests PRIVATE -use_fast_math -O3)
endif()
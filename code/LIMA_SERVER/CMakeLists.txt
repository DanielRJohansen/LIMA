#LIMA_SERVER

find_package(CUDA REQUIRED)

add_executable(mdrun
    mdrun.cpp
)

target_include_directories(mdrun PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(mdrun LIMA_MD)

set_target_properties(mdrun PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(mdrun PROPERTIES CUDA_ARCHITECTURES "89")

if (UNIX)
    target_compile_options(limatests PRIVATE -use_fast_math -O3)
endif()


#set_target_properties(mdrun PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS mdrun DESTINATION .)

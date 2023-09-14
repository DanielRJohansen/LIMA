#LIMA_BASE
find_package(CUDA REQUIRED)

add_library(LIMA_BASE STATIC
    src/Printer.cpp
    src/Bodies.cu
    src/Filehandling.cpp
    src/Forcefield.cu
    src/LimaTypes.cu
    src/Simulation.cu
    src/utilities.cpp
)



target_include_directories(LIMA_BASE PRIVATE ${CUDA_INCLUDE_DIRS})
target_include_directories(LIMA_BASE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

set_target_properties(LIMA_BASE PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_BASE PROPERTIES CUDA_ARCHITECTURES "89")

if (UNIX)
    target_compile_options(LIMA_BASE PRIVATE -use_fast_math -O3)
endif()

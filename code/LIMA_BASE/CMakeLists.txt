#LIMA_BASE 

add_library(LIMA_BASE STATIC
    src/Printer.cpp
    src/Bodies.cu
    src/Filehandling.cpp
    src/Forcefield.cpp
    src/LimaTypes.cpp
    src/Simulation.cpp
    src/utilities.cpp
    src/MoleculeGraph.cpp
    src/MDFiles.cpp
    src/TrrFileInterface.cpp
)

target_include_directories(LIMA_BASE PRIVATE ${CUDA_INCLUDE_DIRS})
target_include_directories(LIMA_BASE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

if (UNIX)
    target_compile_options(LIMA_BASE PRIVATE -use_fast_math -O3)
endif()

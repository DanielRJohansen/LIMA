#LIMA_BASE
cmake_minimum_required(VERSION 3.27 FATAL_ERROR)


add_library(LIMA_BASE STATIC
    src/Printer.cpp
    src/Bodies.cu
    src/Filehandling.cpp
    src/Forcefield.cu
    src/LimaTypes.cu
    src/Simulation.cu
    src/utilities.cpp
)

target_include_directories(LIMA_BASE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


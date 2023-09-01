#LIMA_MD
#cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

find_package(CUDA REQUIRED)


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


target_include_directories(LIMA_MD PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    ${CUDA_INCLUDE_DIRS}
)
target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)


target_link_libraries(LIMA_ENGINE LIMA_BASE)
target_link_libraries(LIMA_ENGINE LIMA_ENGINE)
target_link_libraries(LIMA_ENGINE LIMA_FORCEFIELDMAKER)


set_target_properties(LIMA_MD PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_MD PROPERTIES CUDA_ARCHITECTURES "89")

#LIMA_BASE 

add_library(LIMA_BASE STATIC
    src/Printer.cpp
    src/Bodies.cpp
    src/Filehandling.cpp
    src/Forcefield.cpp
    src/LimaTypes.cpp
    src/Simulation.cpp
    src/utilities.cpp
    src/MoleculeGraph.cpp
    src/MDFiles.cpp
    src/TrrFileInterface.cpp
    src/Statistics.cpp
    src/BoundaryConditionPublic.cpp
    src/MoleculeHull.cpp
    src/TimeIt.cpp
    src/Geometry.cpp
    src/MoleculeUtils.cpp 
    src/Topologyfile.cpp
    src/BoxGrid.cpp
    src/PlotUtils.cpp
)

target_include_directories(LIMA_BASE PRIVATE 
 #   /../dependencies/cereal
)
target_include_directories(LIMA_BASE PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/code/dependencies/cereal/include
    ${CUDAToolkit_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/code/dependencies/
    ${CMAKE_SOURCE_DIR}/code/dependencies/glm/

)

#target_precompile_headers(LIMA_BASE PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:<vector>>" )
if (UNIX)  
    target_link_libraries(LIMA_BASE PUBLIC
        LIMA_DEPENDS
        TBB::tbb
    )
else()
target_link_libraries(LIMA_BASE PUBLIC
        LIMA_DEPENDS
    )
endif()

if (UNIX)
    target_compile_options(LIMA_BASE PRIVATE -use_fast_math -O3)
else()
    target_compile_options(LIMA_BASE PRIVATE)
endif()

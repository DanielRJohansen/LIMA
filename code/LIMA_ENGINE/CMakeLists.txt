#LIMA_ENGINE

# Do NOT add any other .cu files, the engine should be a single compilation unit
# in order to access __constant__ memory
add_library(LIMA_ENGINE STATIC
    src/Engine.cu
    "src/Thermostat.cpp"       # TODO: Make .cpp file
    src/ChargeOcttree.cu
    src/ConvexHullEngine.cu
 "src/SimulationDevice.cpp")

target_include_directories(LIMA_ENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_ENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_ENGINE LIMA_BASE)
#target_link_libraries(LIMA_ENGINE ${PYTHON_LIBRARIES})

if (UNIX)
    target_compile_options(LIMA_ENGINE PRIVATE -use_fast_math -O3)
elseif (WIN32)
       target_compile_options(LIMA_ENGINE PRIVATE 
        $<$<CONFIG:Release>:/use_fast_math>
        $<$<CONFIG:Release>:/O2>
    )
endif()

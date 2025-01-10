#LIMA_CONVEXHULLENGINE

# Do NOT add any other .cu files, the engine should be a single compilation unit
# in order to access __constant__ memory
add_library(LIMA_CONVEXHULLENGINE STATIC
    src/ConvexHullEngine.cu
)

target_include_directories(LIMA_CONVEXHULLENGINE PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_CONVEXHULLENGINE PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_CONVEXHULLENGINE LIMA_BASE LIMA_ENGINE)
#target_link_libraries(LIMA_CONVEXHULLENGINE ${PYTHON_LIBRARIES})

#target_compile_options(LIMA_CONVEXHULLENGINE PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
#                       --generate-line-info
#                       --use_fast_math
#                       --relocatable-device-code=true
#                       >)

if (UNIX)
    target_compile_options(LIMA_CONVEXHULLENGINE PRIVATE -use_fast_math -O3)
elseif (WIN32)
       target_compile_options(LIMA_CONVEXHULLENGINE PRIVATE 
        $<$<CONFIG:Release>:/use_fast_math>
        $<$<CONFIG:Release>:/O2>
    )
endif()

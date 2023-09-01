#LIMA_MD
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


target_include_directories(LIMA_MD PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/../dependencies/GLFW/include
    ${CUDA_INCLUDE_DIRS}    

)
target_include_directories(LIMA_MD PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)



target_link_libraries(LIMA_MD LIMA_BASE)
target_link_libraries(LIMA_MD LIMA_ENGINE)
target_link_libraries(LIMA_MD LIMA_FORCEFIELDMAKER)
target_link_libraries(LIMA_MD ${CMAKE_SOURCE_DIR}/../dependencies/GLFW/lib/glfw3.lib)
#target_include_directories(LIMA_MD PRIVATE ${OPENGL_INCLUDE_DIR})
target_link_libraries(LIMA_MD ${OPENGL_LIBRARIES})


set_target_properties(LIMA_MD PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
set_target_properties(LIMA_MD PROPERTIES CUDA_ARCHITECTURES "89")





# Self-test executable
#add_executable(self_test src/Main.cpp)
#target_link_libraries(self_test LIMA_MD)
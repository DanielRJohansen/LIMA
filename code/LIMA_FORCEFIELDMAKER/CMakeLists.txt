#LIMA_FORCEFIELDMAKER
#cmake_minimum_required(VERSION 3.27 FATAL_ERROR)

#find_package(CUDA REQUIRED)


add_library(LIMA_FORCEFIELDMAKER STATIC
    src/ForcefieldTypes.cpp
    src/LIMA_ForcefieldMaker.cpp
)


#target_include_directories(LIMA_FORCEFIELDMAKER PRIVATE ${CUDA_INCLUDE_DIRS})
target_include_directories(LIMA_FORCEFIELDMAKER PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_ENGINE LIMA_BASE)


#set_target_properties(LIMA_FORCEFIELDMAKER PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
#set_target_properties(LIMA_FORCEFIELDMAKER PROPERTIES CUDA_ARCHITECTURES "89")

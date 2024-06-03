#LIMA_FORCEFIELDMAKER

add_library(LIMA_FORCEFIELDMAKER STATIC    
    "src/Forcefield.cpp"
)


target_include_directories(LIMA_FORCEFIELDMAKER PRIVATE
)
target_include_directories(LIMA_FORCEFIELDMAKER PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(LIMA_FORCEFIELDMAKER LIMA_BASE)


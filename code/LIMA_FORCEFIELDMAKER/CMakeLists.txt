#LIMA_FORCEFIELDMAKER

add_library(LIMA_FORCEFIELDMAKER STATIC
    src/ForcefieldTypes.cpp
    src/LIMA_ForcefieldMaker.cpp
)


target_include_directories(LIMA_FORCEFIELDMAKER PRIVATE
)
target_include_directories(LIMA_FORCEFIELDMAKER PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(LIMA_FORCEFIELDMAKER LIMA_BASE)


# Self-test executable
add_executable(ffm_self_test src/Main.cpp)
target_link_libraries(ffm_self_test LIMA_BASE)
target_link_libraries(ffm_self_test LIMA_FORCEFIELDMAKER)

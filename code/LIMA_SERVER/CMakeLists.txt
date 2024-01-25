#LIMA_SERVER



# First compile MDRUN
add_executable(mdrun
    mdrun.cpp
)

target_include_directories(mdrun PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(mdrun LIMA_MD)

install(TARGETS mdrun DESTINATION .)



# Now compile rest
add_executable(lima
    lima.cpp
)

target_include_directories(lima PRIVATE
    ${CUDA_INCLUDE_DIRS}
)

target_link_libraries(lima LIMA_MD)



#set_target_properties(lima PROPERTIES INSTALL_PREFIX ${CMAKE_SOURCE_DIR})
install(TARGETS lima DESTINATION .)
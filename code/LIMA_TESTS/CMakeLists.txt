#LIMA_TESTS

if(WIN32)
    enable_language(RC)
    set(RESOURCE_FILES
        ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc
    )
endif()

set(SOURCE_FILES
    src/main.cpp
)
add_executable(limatests
    #src/main.cpp
    ${SOURCE_FILES}
    ${RESOURCE_FILES}
)

target_include_directories(limatests PRIVATE
    ${CUDA_INCLUDE_DIRS}
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(limatests LIMA_BASE)
target_link_libraries(limatests LIMA_MD)

if (UNIX)
    target_compile_options(limatests PRIVATE -use_fast_math -O3)
endif()

message(STATUS "\n Resource files: ${RESOURCE_FILES}")
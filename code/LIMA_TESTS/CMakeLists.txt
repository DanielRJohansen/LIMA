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
add_executable(LIMA_TESTS
    ${SOURCE_FILES}
    ${RESOURCE_FILES}
 "src/ForcefieldTests.h")

target_include_directories(LIMA_TESTS PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(LIMA_TESTS LIMA_BASE)
target_link_libraries(LIMA_TESTS LIMA_MD)

if (UNIX)
    target_compile_options(LIMA_TESTS PRIVATE -use_fast_math -O3)
endif()

message(STATUS "\n Resource files: ${RESOURCE_FILES}")
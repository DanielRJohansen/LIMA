#limatest

if(WIN32)
    enable_language(RC)
    set(RESOURCE_FILES
        ${CMAKE_SOURCE_DIR}/resources/logo/logoicon.rc
    )
endif()

set(SOURCE_FILES
    src/main.cpp
)
add_executable(limatest
    ${SOURCE_FILES}
    ${RESOURCE_FILES}
 )

target_include_directories(limatest PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_link_libraries(limatest LIMA_BASE)
target_link_libraries(limatest LIMA_MD)


if (UNIX)
    target_compile_options(limatest PRIVATE -O3)
    install(TARGETS limatest DESTINATION /usr/bin)
elseif (WIN32)
    #target_compile_options(LIMA_TESTS PRIVATE /O2)
endif()

message(STATUS "\n Resource files: ${RESOURCE_FILES}")

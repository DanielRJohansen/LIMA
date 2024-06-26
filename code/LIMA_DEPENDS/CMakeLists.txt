#LIMA_DEPENDS
# Manually compiled dependencies for LIMA

add_library(LIMA_DEPENDS STATIC
    glew-2.1.0/src/glew.c
)

target_include_directories(LIMA_DEPENDS PUBLIC
    glew-2.1.0/include/
)

target_compile_definitions(LIMA_DEPENDS PUBLIC GLEW_STATIC)

set_target_properties(LIMA_DEPENDS PROPERTIES LINKER_LANGUAGE CXX)
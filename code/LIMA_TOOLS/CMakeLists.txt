#LIMA_TOOLS

add_library(LIMA_TOOLS STATIC
    src/SimulationBuilder.cpp
    src/CompoundBuilder.cpp
    src/BondgroupFactory.cpp
)

target_include_directories(LIMA_TOOLS PRIVATE
${CMAKE_CURRENT_SOURCE_DIR}/src
)

target_include_directories(LIMA_TOOLS PUBLIC
${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(LIMA_TOOLS LIMA_BASE)
target_link_libraries(LIMA_TOOLS LIMA_MD)

if (UNIX)
    target_compile_options(LIMA_TOOLS PRIVATE -use_fast_math -O3)
endif()

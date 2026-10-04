#LIMA_DEPENDS
# Manually compiled dependencies for LIMA

add_library(LIMA_DEPENDS STATIC
    glew-2.1.0/src/glew.c

    qhull/src/libqhullcpp/Coordinates.cpp
    qhull/src/libqhullcpp/PointCoordinates.cpp
    qhull/src/libqhullcpp/Qhull.cpp
    qhull/src/libqhullcpp/QhullFacet.cpp
    qhull/src/libqhullcpp/QhullFacetList.cpp
    qhull/src/libqhullcpp/QhullFacetSet.cpp
    qhull/src/libqhullcpp/QhullHyperplane.cpp
    qhull/src/libqhullcpp/QhullPoint.cpp
    qhull/src/libqhullcpp/QhullPoints.cpp
    qhull/src/libqhullcpp/QhullPointSet.cpp
    qhull/src/libqhullcpp/QhullQh.cpp
    qhull/src/libqhullcpp/QhullRidge.cpp
    qhull/src/libqhullcpp/QhullSet.cpp
    qhull/src/libqhullcpp/QhullStat.cpp
    qhull/src/libqhullcpp/QhullUser.cpp
    qhull/src/libqhullcpp/QhullVertex.cpp
    qhull/src/libqhullcpp/QhullVertexSet.cpp
    qhull/src/libqhullcpp/RboxPoints.cpp
    qhull/src/libqhullcpp/RoadError.cpp
    qhull/src/libqhullcpp/RoadLogEvent.cpp

    qhull/src/libqhull_r/qset_r.c
    qhull/src/libqhull_r/poly_r.c
    qhull/src/libqhull_r/poly2_r.c
    qhull/src/libqhull_r/merge_r.c
    qhull/src/libqhull_r/geom2_r.c
    qhull/src/libqhull_r/geom_r.c
    qhull/src/libqhull_r/stat_r.c
    qhull/src/libqhull_r/mem_r.c
    qhull/src/libqhull_r/random_r.c
    qhull/src/libqhull_r/usermem_r.c
    qhull/src/libqhull_r/io_r.c
    qhull/src/libqhull_r/libqhull_r.c
    qhull/src/libqhull_r/user_r.c
    qhull/src/libqhull_r/rboxlib_r.c
    qhull/src/libqhull_r/global_r.c

    mdtraj/xtc/src/xdrfile.c
    mdtraj/xtc/src/xdrfile_trr.c
    mdtraj/xtc/src/xdrfile_xtc.c
    mdtraj/xtc/src/xdr_seek.c

    imgui/imgui.cpp
    imgui/imgui_draw.cpp
    imgui/imgui_tables.cpp
    imgui/imgui_widgets.cpp
    imgui/imgui_demo.cpp
    imgui/backends/imgui_impl_glfw.cpp
    imgui/backends/imgui_impl_opengl3.cpp
)

target_include_directories(LIMA_DEPENDS PUBLIC
    ${CMAKE_SOURCE_DIR}/code/LIMA_DEPENDS/imgui/
    ${CMAKE_SOURCE_DIR}/code/LIMA_DEPENDS/imgui/backends
    glew-2.1.0/include/

    qhull/src

    mdtraj/xtc/include/
    ${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/include

)

if (NOT UNIX)
    set(GLFW_LIBRARIES "${CMAKE_SOURCE_DIR}/code/dependencies/GLFW/lib/glfw3.lib")
else()
    set(GLFW_LIBRARIES glfw)
endif()
target_link_libraries(LIMA_DEPENDS PUBLIC
    ${GLFW_LIBRARIES}
    ${OPENGL_LIBRARIES}
)

# GLEW_NO_GLU: LIMA does not use GLU, so dont require its headers or library
target_compile_definitions(LIMA_DEPENDS PUBLIC GLEW_STATIC GLEW_NO_GLU)

set_target_properties(LIMA_DEPENDS PROPERTIES LINKER_LANGUAGE CXX)
# Unity build in hand-picked groups (see LIMA_UNITY): qhull's many small sources share headers and are cheap to merge,
# split in two groups to compile in parallel. imgui's sources are each large, so they stay separate, as do glew.c and
# xdrfile (whose sources redefine each other's types)
if (LIMA_UNITY)
    set_target_properties(LIMA_DEPENDS PROPERTIES UNITY_BUILD ON UNITY_BUILD_MODE GROUP)
    get_target_property(dependsSources LIMA_DEPENDS SOURCES)
    foreach(source ${dependsSources})
        if (source MATCHES "libqhull_r/(qset|poly|poly2|merge|geom2)_r")
            set_source_files_properties(${source} PROPERTIES UNITY_GROUP qhull_r_a)
        elseif (source MATCHES "libqhull_r/")
            set_source_files_properties(${source} PROPERTIES UNITY_GROUP qhull_r_b)
        elseif (source MATCHES "libqhullcpp/")
            set_source_files_properties(${source} PROPERTIES UNITY_GROUP qhullcpp)
        endif()
    endforeach()
endif()

# LIMA does not use the imgui demo window, this compiles it to empty stubs
set_source_files_properties(imgui/imgui_demo.cpp PROPERTIES COMPILE_DEFINITIONS IMGUI_DISABLE_DEMO_WINDOWS)

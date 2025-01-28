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
)

target_include_directories(LIMA_DEPENDS PUBLIC
    glew-2.1.0/include/

    qhull/src

    mdtraj/xtc/include/
)


target_compile_definitions(LIMA_DEPENDS PUBLIC GLEW_STATIC)

set_target_properties(LIMA_DEPENDS PROPERTIES LINKER_LANGUAGE CXX)
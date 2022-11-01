
# Build our own version of the library.
SET(VOROPLUSPLUS_SOURCE_FILES
    src/voro++/cell.cc 
    src/voro++/common.cc 
    src/voro++/container.cc
    src/voro++/unitcell.cc
    src/voro++/v_compute.cc
    src/voro++/c_loops.cc
    src/voro++/v_base.cc
    src/voro++/wall.cc
    src/voro++/pre_container.cc
    src/voro++/container_prd.cc
    )

# Build library.
ADD_LIBRARY(VoroPlusPlus STATIC ${VOROPLUSPLUS_SOURCE_FILES}) 

TARGET_INCLUDE_DIRECTORIES(VoroPlusPlus INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")

INSTALL(TARGETS VoroPlusPlus DESTINATION lib                         
        PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_WRITE GROUP_EXECUTE WORLD_READ WORLD_WRITE WORLD_EXECUTE)



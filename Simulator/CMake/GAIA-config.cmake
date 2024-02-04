cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

message( "Adding GAIA." )


find_package(Eigen3 REQUIRED)
find_package(MeshFrame2 REQUIRED PATHS ${CMAKE_CURRENT_LIST_DIR}/../3rdParty/MeshFrame2/MeshFrame/cmake)
find_package(embree 3.0 REQUIRED)
find_package(CuMatrix REQUIRED PATHS ${CMAKE_CURRENT_LIST_DIR}/../3rdParty/CuMatrix/cmake)

add_subdirectory ("../3rdParty/cmake-git-version-tracking" ${CMAKE_CURRENT_BINARY_DIR}/cmake-git-version-tracking)

set (GAIA_ROOT ${CMAKE_CURRENT_LIST_DIR}/..)

	   
option (BUILD_PBD
       "Build PBD volumetric simulation modules." OFF)

	   

set(THIRD_PARTY_INCLUDE_DIRS
        ${EIGEN3_INCLUDE_DIR}
		${MESHFRAME_INCLUDE_DIR}
		${EMBREE_INCLUDE_DIRS}
		${CU_MATRIX_INCLUDE_DIR}
		"../3rdParty/oneTBB/include"
        )


set(GAIA_INCLUDE_DIRS
	"../Modules"
	${THIRD_PARTY_INCLUDE_DIRS}
	)
	

file(GLOB GAIA_PBD_SRCS
	"../Modules/PBD/*.h"
	"../Modules/PBD/*.cpp"
	"../Modules/PBD/*.cu"
	"../Modules/PBD/*.cuh"
)

file(GLOB GAIA_SRCS
	"../Modules/TetMesh/**.h"
	"../Modules/TetMesh/**.cpp"
	"../Modules/TriMesh/**.h"
	"../Modules/TriMesh/**.cpp"
	"../Modules/Parallelization/*.h"
	"../Modules/Parallelization/*.cpp"
	"../Modules/Parameters/*.h"
	"../Modules/Parameters/*.cpp"
	"../Modules/Parser/*.h"
	"../Modules/Parser/*.cpp"
	"../Modules/VersionTracker/*.h"
	"../Modules/VersionTracker/*.cpp"
	"../Modules/IO/*.h"
	"../Modules/IO/*.cpp"
	"../Modules/CollisionDetector/*.h"
	"../Modules/CollisionDetector/*.cpp"
	"../Modules/SpatialQuery/*.h"
	"../Modules/SpatialQuery/*.cpp"
	"../Modules/Framework/*.h"
	"../Modules/Framework/*.cpp"
	"../Modules/common/math/constants.cpp"

	${MESHFRAME_SOURCE_CPP_UTILITY}
)

list(REMOVE_ITEM GAIA_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/../Modules/CollisionDetector/TetMeshContactDetector.h" "${CMAKE_CURRENT_SOURCE_DIR}/../Modules/CollisionDetector/TetMeshContactDetector.cpp")


if (BUILD_PBD)
message("GAIA: Build with PBD volumetric object simulation components!\n")
SET (GAIA_SRCS 
	${GAIA_SRCS}
	${GAIA_PBD_SRCS}
)
endif (BUILD_PBD)

set (GAIA_LIBRARY
	${CU_MATRIX_LIBS}
	${EMBREE_LIBRARY}
	${embree_DIR}/../../tbb12.lib
	cmake_git_version_tracking
	)

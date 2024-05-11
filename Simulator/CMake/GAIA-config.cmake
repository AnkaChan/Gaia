cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

message( "Adding GAIA." )


find_package(Eigen3 REQUIRED)
find_package(MeshFrame2 REQUIRED PATHS ${CMAKE_CURRENT_LIST_DIR}/../3rdParty/MeshFrame2/MeshFrame/cmake)
find_package(embree 3.0 REQUIRED)
find_package(CuMatrix REQUIRED PATHS ${CMAKE_CURRENT_LIST_DIR}/../3rdParty/CuMatrix/cmake)
find_package(TBB REQUIRED)
set (GAIA_ROOT ${CMAKE_CURRENT_LIST_DIR}/..)
	
add_subdirectory ("${CMAKE_CURRENT_LIST_DIR}/../3rdParty/cmake-git-version-tracking" ${CMAKE_CURRENT_BINARY_DIR}/cmake-git-version-tracking)


option (BUILD_VBD
       "Build VBD modules." ON)
	   
option (BUILD_VBD_Cloth
       "Build VBD Cloth modules." ON)
	   
option (BUILD_PBD
       "Build PBD volumetric simulation modules." OFF)
	   
option (BUILD_GUI
       "Build Gui." ON)
	   
option (BUILD_Collision_Detector
       "Build Collision Detectors Modules." ON)
	   
set(GAIA_DEFINITIONS)



set(THIRD_PARTY_INCLUDE_DIRS
        ${EIGEN3_INCLUDE_DIR}
		${MESHFRAME_INCLUDE_DIR}
		${CU_MATRIX_INCLUDE_DIR}
        )
		
if (BUILD_GUI)
	add_subdirectory ("${CMAKE_CURRENT_LIST_DIR}/../3rdParty/polyscope" ${CMAKE_CURRENT_BINARY_DIR}/polyscope)
	# set(THIRD_PARTY_INCLUDE_DIRS
		# ${THIRD_PARTY_INCLUDE_DIRS}
		# "${CMAKE_CURRENT_LIST_DIR}/../3rdParty/polyscope/include"
	# )
else()
	set(GAIA_DEFINITIONS
		${GAIA_DEFINITIONS}
		GAIA_NO_GUI
	)
endif (BUILD_GUI)

if (BUILD_Collision_Detector)
	message("GAIA: Build with Collision Detectors!\n")
	SET (THIRD_PARTY_INCLUDE_DIRS 
		${THIRD_PARTY_INCLUDE_DIRS}
		${EMBREE_INCLUDE_DIRS}
	)
endif (BUILD_Collision_Detector)


set(GAIA_INCLUDE_DIRS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules"
	${THIRD_PARTY_INCLUDE_DIRS}
	)
	
file(GLOB GAIA_PBD_SRCS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/PBD/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/PBD/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/PBD/*.cu"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/PBD/*.cuh"
)

file(GLOB GAIA_VBD_SRCS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VBD/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VBD/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VBD/*.cu"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VBD/*.cuh"
)


file (GLOB GAIA_COLLISION_SRCS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/BVH/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/BVH/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/BVH/*.cu"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/SpatialQuery/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/SpatialQuery/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/SpatialQuery/*.cu"
)

file(GLOB GAIA_SRCS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/TetMesh/**.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/TetMesh/**.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parallelization/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parallelization/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/SolverUtils/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/SolverUtils/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parameters/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parameters/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parser/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Parser/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VersionTracker/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/VersionTracker/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/IO/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/IO/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Framework/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Framework/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/common/math/constants.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/Viewer/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/CollisionDetector/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/CollisionDetector/*.cpp"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/TriMesh/**.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/TriMesh/**.cpp"
	${MESHFRAME_SOURCE_CPP_UTILITY}
)

file(GLOB GAIA_CLOTH_SRCS

)

file(GLOB GAIA_COLORING_SRCS
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/GraphColoring/*.h"
	"${CMAKE_CURRENT_LIST_DIR}/../Modules/GraphColoring/*.cpp"
)

list(REMOVE_ITEM GAIA_COLLISION_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_CURRENT_LIST_DIR}/../Modules/CollisionDetector/TetMeshContactDetector.h" "${CMAKE_CURRENT_SOURCE_DIR}/${CMAKE_CURRENT_LIST_DIR}/../Modules/CollisionDetector/TetMeshContactDetector.cpp")

if (BUILD_VBD)
message("GAIA: Build with VBD components!\n")
SET (GAIA_SRCS 
	${GAIA_SRCS}
	${GAIA_VBD_SRCS}
)
endif (BUILD_VBD)

if (BUILD_PBD)
message("GAIA: Build with PBD volumetric object simulation components!\n")
SET (GAIA_SRCS 
	${GAIA_SRCS}
	${GAIA_PBD_SRCS}
)
endif (BUILD_PBD)


if (BUILD_Collision_Detector)
message("GAIA: Build with Collision Detector components!\n")
SET (GAIA_SRCS 
	${GAIA_SRCS}
	${GAIA_COLLISION_SRCS}
)
endif (BUILD_Collision_Detector)

if (BUILD_VBD_Cloth)
message("GAIA: Build with VBD Cloth simulatio components!\n")
SET (GAIA_SRCS 
	${GAIA_SRCS}
	${GAIA_CLOTH_SRCS}
)
endif (BUILD_VBD_Cloth)



set (GAIA_LIBRARY
	${CU_MATRIX_LIBS}
	${EMBREE_LIBRARY}
	cmake_git_version_tracking
	TBB::tbb
	)

if (BUILD_GUI)
	set (GAIA_LIBRARY
		${GAIA_LIBRARY}
		polyscope
	)
endif(BUILD_GUI)
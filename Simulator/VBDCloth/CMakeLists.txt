cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
# add_compile_options($<$<CXX_COMPILER_ID:MSVC>:/MP>)

project(GAIA_VBDCloth LANGUAGES CXX CUDA)

## Use C++11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)   

set(BUILD_PBD OFF)
set(BUILD_VBD OFF)
set(BUILD_VBD_Cloth ON)
find_package(GAIA REQUIRED)

include_directories(
	${GAIA_INCLUDE_DIRS}
	${CMAKE_CURRENT_LIST_DIR}
)

file(GLOB SRC
    "*.h"
    "*.cpp"
	"*.c"
	)

add_executable(GAIA_VBDCloth 
	${SRC}
	${GAIA_SRCS}
	${CNPY_SRC}
)

target_include_directories(GAIA_VBDCloth PUBLIC ${CMAKE_CURRENT_LIST_DIR}/cnpy)

option(USE_MKL "whether to use intel mkl to solve linear systems, if not use eigen" FALSE)
option(USE_DOUBLE "whether to use double to solve linear system" OFF)


if(USE_MKL)
	set(MKL_INTERFACE "lp64" CACHE STRING "" FORCE)
	set(MKL_INTERFACE_FULL "intel_lp64" CACHE STRING "" FORCE)
	find_package(MKL CONFIG REQUIRED)
	target_compile_options(GAIA_VBDCloth PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_compile_definitions(GAIA_VBDCloth PUBLIC EIGEN_USE_MKL_ALL USE_MKL)
	target_link_libraries(GAIA_VBDCloth MKL::MKL)
	target_include_directories(GAIA_VBDCloth PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
endif()


target_compile_definitions(GAIA_VBDCloth PUBLIC ${GAIA_DEFINITIONS})

target_compile_options(GAIA_VBDCloth PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
                       --extended-lambda
					   --default-stream per-thread
                       >)

target_link_libraries(GAIA_VBDCloth ${GAIA_LIBRARY})
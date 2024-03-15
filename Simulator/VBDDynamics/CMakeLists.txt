cmake_minimum_required(VERSION 3.13 FATAL_ERROR)

find_package(CUDAToolkit 11 REQUIRED)

project(VBDDynamics LANGUAGES CXX CUDA)

## Use C++11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)   

set(BUILD_EBD OFF)
set(BUILD_VBD ON)
set(BUILD_PBD OFF)
set(BUILD_APAP_Cloth OFF)
set(BUILD_VBD_Cloth OFF)
include(../cmake/GAIA-config.cmake)

include_directories(
	${GAIA_INCLUDE_DIRS}
	
)

file(GLOB SRC
    "*.h"
    "*.cpp"
	"*.c"
	"*.cu"
	)

add_executable(VBDDynamics 
	${SRC}
	${GAIA_SRCS}
)


target_compile_options(VBDDynamics PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
                       --extended-lambda
					   --default-stream per-thread
                       >)
#target_compile_definitions(P10_VBDDynamics PUBLIC USE_IPC_FRICTION USE_IPC_BOUNDARY_FRICTION)

target_link_libraries(VBDDynamics ${GAIA_LIBRARY})

option(USE_MKL "whether to use intel mkl to solve linear systems, if not use eigen" OFF)
option(USE_DOUBLE "whether to use double to solve linear system" OFF)
option(PSD_FILTERING "whether to use psd projection" ON)

if(USE_MKL)
	set(MKL_INTERFACE "lp64" CACHE STRING "" FORCE)
	set(MKL_INTERFACE_FULL "intel_lp64" CACHE STRING "" FORCE)
	find_package(MKL CONFIG REQUIRED)
	target_compile_options(VBDDynamics PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
	target_compile_definitions(VBDDynamics PUBLIC EIGEN_USE_MKL_ALL USE_MKL)
	target_link_libraries(VBDDynamics MKL::MKL)
	target_include_directories(VBDDynamics PUBLIC $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
endif()

if(USE_DOUBLE)
	target_compile_definitions(VBDDynamics PUBLIC USE_DOUBLE)
endif()

if(PSD_FILTERING)
	target_compile_definitions(VBDDynamics PUBLIC PSD_FILTERING)
endif()

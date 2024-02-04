cmake_minimum_required(VERSION 3.5 FATAL_ERROR)
project(PBDDynamics LANGUAGES CXX CUDA)

## Use C++11
set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CUDA_STANDARD 11)   

set(BUILD_PBD ON)
include(../cmake/GAIA-config.cmake)

include_directories(
	${GAIA_INCLUDE_DIRS}
	
)

file(GLOB SRC
    "*.h"
    "*.cpp"
	"*.c"
	)

add_executable(PBDDynamics 
	${SRC}
	${GAIA_SRCS}
)


target_compile_options(PBDDynamics PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:
                       --extended-lambda
					   --default-stream per-thread
                       >)

target_link_libraries(PBDDynamics ${GAIA_LIBRARY})
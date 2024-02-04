# Gaia

<a href=""><img src="https://github.com/AnkaChan/Gaia/tree/main/teaser.gif" height="192px"></a>

The Gaia engine is a C++ codebase primarily designed for physics-based simulations. It can be compiled as a standalone simulator or integrated into other applications as a third-party module.
It provides a set of useful tools, including a powerful triangular/tet mesh data structure, a convenient parameter IO module, a set of efficient collision detectors, and a virtual physics framework that can be extended to support all sorts of solvers.


Gaia is engineered to enhance the efficiency of both developers and the hardware it operates on. Occasionally, compromises are necessary since optimizing for one can adversely affect the other. It may not have the smartest design in every way. However, it's assured that Gaia is free from any stupid designs.

## Installation

### Get the Code
Please download the code using  the following command:
```
git clone git@github.com:AnkaChan/Gaia.git --recursive
```

### Dependencies
This Algorithm has the following dependencies:
- [MeshFrame2](https://github.com/AnkaChan/MeshFrame2): geometric core library for mesh processing (already included as submodule)
- [CuMatrix](https://github.com/AnkaChan/CuMatrix/tree/main): for geometry and matrix computation (already included as submodule)
- [cmake-git-version-tracking](https://github.com/andrew-hardin/cmake-git-version-tracking): for tracking the git version infos when running experiments (already included as a submodule)
- [OneTBB](https://github.com/oneapi-src/oneTBB) (already included as a submodule)
- Eigen3 (tested with 3.4.0)
- Embree (tested with 3.13.1) Known issues exist for version >= 3.13.4.

Ensure you install the specified versions of Eigen3 and Embree, and then add the environment variables "Eigen3_DIR" and "embree_DIR", setting their values to the respective config.cmake paths. This step is necessary for CMake to successfully locate them.

### Use GAIA as a Standalone Simulator
Currently, a XPBD based simulator is provided. The CMakelists.txt of the standalone simulator is located at Gaia/Simulator/PBDDynamics. To use the simulator, you need to
use CMake to construct and build it. As of now, compatibility and testing have been conducted exclusively on Windows with Visual Studio.

### Use GAIA as a Module
Gaia has its own GAIA-config.cmake located at Gaia/Simulator/CMake. To incorperate Gaia as a module, you need to add ```include([PATH_TO_GAIA]/Simulator/cmake/GAIA-config.cmake)``` to your project's CMakeLists, then add those commands to link Gaia to your project:
```
include_directories(
	${GAIA_INCLUDE_DIRS}
)
SET (YOUR_SRCS 
	${YOUR_SRCS}
	${GAIA_SRCS}
)
add_executable(YouApplication 
	${YOUR_SRCS}
)
target_link_libraries(YouApplication ${GAIA_LIBRARY})
```
An alternative way to make GAIA visible to CMake is defining an Environment Variable named "GAIA_DIR", setting its value to [PATH_TO_GAIA]/Simulator/cmake/, and use CMake command ```find_package(GAIA)```.

The Gaia engine provides various modules that can be enabled or disabled, applicable in both Standalone and Module modes. This functionality is achieved through a series of CMake options:
```set(GAIA_OPTION ON/OFF)```.
You need to put this option before the command ```include([PATH_TO_GAIA]/Simulator/cmake/GAIA-config.cmake)``` or ```find_package(GAIA)``` to make it effective.
Those options determine the source files included into ```${GAIA_SRCS}```, as well as the static libraries included into ```${GAIA_LIBRARY}```.
Currently, only one option is avaliable, that is ```BUILD_PBD```.

## Run Gaia

GAIA simulators accept three positional command-line arguments along with several optional keyword arguments. Here's the basic syntax for running a GAIA simulator:
```
Path-to-PBDDynamics.exe Models.json Parameters.json output-folder -R [PATH-to-Gaia-Repository]
```

The first argument, "Models.json", specifies the details of the models to be simulated, while "Parameters.json" contains the physics parameters for the simulation. The third argument, "output-folder", designates the directory where the simulation results will be stored. The keyword argument "-R [PATH-to-Gaia-Repository]" is used to replace the placeholder "${REPO_ROOT}" found in "Models.json" and "Parameters.json" with the actual path to the Gaia repository, simplifying the loading process. To discover more command-line options, execute:
```PBDDyanmics -h```.

Given the complexity of manually creating "Models.json" and "Parameters.json", GAIA offers some Python scripts that can automatically generate the simulation command for you. Those scripts covers most of the experiments included in my paper: "Shortest Path to Boundary for Self Intersecting Meshes".

## Common Problems

1. "tbb12.lib" is not found.  
The reason is that Embree has already had a compilation of tbb, named tbb.lib. However, OneTBB is asking for a file called tbb12.lib. An easy solution is to duplicate that tbb.lib from embree, name it "tbb12.lib" and put it in the same repository.
2. "Host key verification failed. fatal: Could not read from remote repository."  
This issue could be that Github isn't present in your ~/.ssh/known_hosts file.
Append GitHub to the list of authorized hosts:
```ssh-keyscan -H github.com >> ~/.ssh/known_hosts```


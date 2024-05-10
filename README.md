# Gaia

[\[Discord\]](https://discord.gg/BRRaSmRpkm "Gaia Engine")

<img alt="teaser.gif" src="https://github.com/AnkaChan/Gaia/blob/main/teaser.gif?raw=true" data-hpc="true" class="Box-sc-g0xbh4-0 kzRgrI" height="256px">

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
- [OneTBB](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#onetbb) (Need to be installed, tested with 2021.12.0)
- Eigen3 (Need to be installed, tested with 3.4.0)
- Embree (Need to be installed, tested with 3.13.1, not compatible with Embree 4)'
- [polyscope](https://github.com/nmwsharp/polyscope): needed if the BUILD_GUI option for CMAKE is on. (already included as submodule)

Make you installed OneTBB, Eigen3 and Embree. Then add the environment variables "Eigen3_DIR" and "embree_DIR", setting their values to the respective config.cmake paths. This step is necessary for CMake to successfully locate them.

### Use GAIA as a Standalone Simulator
Currently, both VBD (Vertex Block Descent) and XPBD (Extended Position Based Dynamics) based simulatora re provided. The CMakelists.txt of the standalone simulators are located at Gaia/Simulator/VBDDynamics and Gaia/Simulator/PBDDynamics, respectively. To use the simulator, you need to
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

## Gaia's VBD (Vertex Block Descent) Simulator  

GAIA simulators accept three positional command-line arguments along with several optional keyword arguments. Here's the basic syntax for running a GAIA simulator:
```
Path-to-VBDDynamics.exe Models.json Parameters.json output-folder -R [PATH-to-Gaia-Repository]
```
I have already prepared some parameters, located at /Simulator/VBDDynamics/ParameterGen/Parameters/. Run the following command and you should see a simulation of 32 models dropping into a box:
```
Path-to-VBDDynamics.exe [PATH-to-Gaia-Repository]/Simulator/VBDDynamics/ParameterGen/Parameters/S01_Experiment_HybridModelsDrop_sequentiallyAppea/S01_Experiment_HybridModelsDrop_sequentiallyAppear/Models.json [PATH-to-Gaia-Repository]/Simulator/VBDDynamics/ParameterGen/Parameters/S01_Experiment_HybridModelsDrop_sequentiallyAppea/S01_Experiment_HybridModelsDrop_sequentiallyAppear/Parameters.json noOutput -R [PATH-to-Gaia-Repository] --gui
```

Please remember to replace "[PATH-to-Gaia-Repository]" with the absolute path of Gaia repository.

The first argument, namely the file "Models.json", specifies the details of the models to be simulated, while "Parameters.json" contains the physics parameters for the simulation. The third argument, "output-folder", designates the directory where the simulation results will be stored. If it is set to "noOuput", the simulator will not save the simulation results.
The keyword argument "-R [PATH-to-Gaia-Repository]" is used to replace the placeholder "${REPO_ROOT}" found in "Models.json" and "Parameters.json" with the actual path to the Gaia repository, simplifying the loading process. To discover more command-line options, execute:
```VBDDyanmics -h```.

Given the complexity of manually creating "Models.json" and "Parameters.json", GAIA offers the Python scripts that can generates those simulation parameters. Those scripts covers most of the experiments included in my paper: "Vertex Block Descent". Please see: /Simulator/VBDDynamics/ParameterGen/ for those commands.

Additionally, you can run the experiments using these Python scripts by setting the "run" variable to True. This offers a more intuitive way to tune parameters and automate your experiments.  To ensure the Python script recognizes your environment, you need to add the appropriate attribute to the "machines" dictionary in /Simulator/VBDDynamics/ParameterGen/M02_GenRunningParameters.py:
```
    "EnvironmentName": {
        "binaryFile": Absolute-Path-to-VBDDynamics.exe,
        "RepoPath": str(pathlib.Path(__file__).parent.parent.parent.parent.resolve()), #don't change this
        "OutputDirs": [
            Your-Preferred-Output-Path,
        ]
    },
```
Then set the "machineName" variable to your environment's name. The Python script will be able to execute the simulation.

## Gaia's PBD Simulator  
The PBD simulator uses commands similar to those of the VBD simulator. However, please do not input VBD's parameters into the PBD simulator and vice versa. Instead, use the parameters located at \Simulator\PBDDynamics\ParameterGen. The syntax for running the Gaia PBD simulator is similar to that of the VBD simulator:
```
Path-to-VBDDynamics.exe Models.json Parameters.json output-folder -R [PATH-to-Gaia-Repository]
```
## Tetmesh and Trimesh Coloring
A coloring toolkit is located at:Gaia\Simulator\GraphColoring.
Exemplar usage:
```
GraphColoring.exe model.t model.vertexColoring.json -v -b
```
"-v" means it's a vertex coloring. "-b" means it will try to balance the number of vertices in each color.

## Common Problems
1. "Host key verification failed. fatal: Could not read from remote repository."  
This issue could be that Github isn't present in your ~/.ssh/known_hosts file.
Append GitHub to the list of authorized hosts:
```ssh-keyscan -H github.com >> ~/.ssh/known_hosts```


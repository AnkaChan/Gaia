# Gaia
Gaia Physics Engine

Gaia engine is a C++ code base designed mainly for physics based simulation. It provides a set of useful tools, including a powerful triangular/tet mesh data structure, a convenient parameter IO module, a set of efficient collision detectors, and a virtual physics framework that can be extended to support all sorts of solvers.

Gaia is engineered to enhance the efficiency of both developers and the hardware it operates on. Occasionally, compromises are necessary since optimizing for one can adversely affect the other. It may not have the smartest design in every way. However, it's assured that Gaia is free from any stupid designs.

## Installation

### Get the Code
Please download the code using  the following command:
```
git clone git@github.com:AnkaChan/Gaia.git --recursive
```

### Dependencies
This Algorithm has the following dependencies:
- [MeshFrame2](https://github.com/AnkaChan/MeshFrame2): for mesh processing (included as submodule)
- [CuMatrix](https://github.com/AnkaChan/CuMatrix/tree/main): for geometry and matrix computation (included as submodule)
- [cmake-git-version-tracking](https://github.com/andrew-hardin/cmake-git-version-tracking): for geometry and matrix computation (included as submodule)
- OneTBB (included as a submodule)
- Eigen3 (tested with 3.4.0)
- Embree (tested with 3.13.1) Known issues exist for version >= 3.13.4.

You need to install Eigen3 and Embree with the required version and add attributes: "Eigen3_DIR" and "embree_DIR" whose values are their corresponding config.cmake path to your environment variable to allow CMake to find them.

### Compile
Use CMake to build the project and compile it. Currently, the code is only tested with Windows & Visual Studio. 

### Common Problems

1. "tbb12.lib" is not found.  
The reason of this bug is that Embree has already had a compilation of tbb, named tbb.lib. However, OneTBB is asking for a file called tbb12.lib. An easy solution is to duplicate that tbb.lib from embree, name it "tbb12.lib" and put it in the same repository.
2. "Host key verification failed. fatal: Could not read from remote repository."  
This issue could be that Github isn't present in your ~/.ssh/known_hosts file.
Append GitHub to the list of authorized hosts:
```ssh-keyscan -H github.com >> ~/.ssh/known_hosts```


## Test and Run
After compilation, it should give you a binary file called: Shortest-Path-Test.
Run it using the following command:

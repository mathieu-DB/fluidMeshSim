# Fluid Simulation on a Mesh

##Build
C++ and cmake are necessary
Clone this repository 
```
git clone https://github.com/mathieu-DB/fluidMeshSim.git
```
libigl is not included in the repo
```
cd fluidMeshSim
git clone https://github.com/libigl/libigl.git
```
Create necessary files
```
mkdir build
cmake ..
```
If your system created a Visual Studio project, you can
* Open the project and BUILD (not run) from VS
or
* Configure msbuild in your path and run 
```
msbuild fluidmesh.sln
```
If your system created another type of build (makefile, ninja.build, etc...), use the proper command to build the project.

Locate your executable file and run the project with the provided `camel_b.obj` file
```
fluidMeshSim.exe [PATH TO CAMEL_B.OBJ]
```

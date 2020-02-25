## Hilltests

This folder contains the files to run the hilltest for configuration ii. 

The standard setup is a Crank-Nicholson time discretization and SUPG stabilisation in space. 
However, if you want to test the other spatial stabilisation schemes (CIP, BR, SOLD, DG) checkout the different moving_surface_Stabilisationabbreviation.py files.
Here, the msh file is directly updated when a new 3D mesh is generated.
In this example code, the mesh has 4 vertical layers and by using the msh approach, this number is fixed. 
Find all this in the folder: Standard_Setup

If you want to use adaptive mesh refinement, test the version in folder: Adaptive_Mesh_Refinement.
Here, a new stl file is created with different spatial resolution according to the refinement algorithm. 
Hence, the number of vertical layers in not fixed and also the mesh is not structured in the vertical.

A 2nd order Runge Kutta time discretization option is provided in folder: Runge_Kutta.
Here, the msh file is directly updated when a new 3D mesh is generated.
In this example code, the mesh has 4 vertical layers and by using the msh approach, this number is fixed.


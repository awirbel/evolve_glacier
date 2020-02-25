# evolve_glacier

This repository contains a development version of a FEM model for the evolution of the free surface of  a glacier due to ice flow and mass balance rate fully accounting for the constraint of **S >= B**. (**S** free surface elevation and **B** bed elevation)

The model consists of different modules that automatizes ice flow and free surface evolution computations and generation of updated computational meshes. 

Ice flow computations are performed with the FEM ice flow model [icetools](https://github.com/alexjarosch/icetools).

A manuscript describing the underlying phyics, numerical details and glaciological relevance is submitted to GMD (Geoscientific Model Development). 


To run the model, the simplest way is to run a singularity container that includes a working installation of FEniCS v2016.2 and * [Gmsh](http://gmsh.info/).

## Requirements
* [singularity v2.6](https://www.sylabs.io/guides/2.6/user-guide/installation.html)

A singularity container with FEniCS v2016.2 and gmsh can be found as an asset with the first release of this code. 

To run the model execute the following line in a shell:
```shell
singularity exec pathtocontainer/fenics_v2.simg ./run_test01.sh > log01.txt
```

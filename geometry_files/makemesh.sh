#!/bin/bash

# create box STL file
singularity exec pathtocontainer/fenics_v2.simg gmsh-aj -2 -format stl surface_shell_domain.geo

# create initial 3D geometry stl file from surf and bed tifs and box STL file
python process_vol_tif_ail.py

# generate 3D msh file
singularity exec pathtocontainer/fenics_v2.simg gmsh-aj -3 glacier_0000.geo

# convert to xml file
singularity exec pathtocontainer/fenics_v2.simg dolfin-convert glacier_0000.msh glacier_0000.xml

# create computational 2D surf_file
singularity exec pathtocontainer/fenics_v2.simg gmsh-aj -2 surface_2D.geo
singularity exec pathtocontainer/fenics_v2.simg dolfin-convert surface_2D.msh surface_2D.xml

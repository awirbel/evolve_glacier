#!/bin/bash
set -e

# Define data directory
name="glacier_"
name2=".xml"
name_dir="/workingdirectory/data/"

res3=$(date +%s.%N)

# Load configuration and mesh
python cfg.py
python read_mesh.py

# Compute 3D velocity field for given mesh
res1=$(date +%s.%N)
mpirun -np 4 python ice_dynamics3D.py
res2=$(date +%s.%N)
awk "BEGIN {printf \"Iceflow runtime [seconds]: %.5f\n\", $res2-$res1}"

# Compute free surface evolution
res1=$(date +%s.%N)
cp surface_start1_shell.geo $name_dir$"surface_shell_domain.geo"
cp surface_start1_2D.geo $name_dir$"surface_2D_0001.geo"
python moving_surface_ini.py
res2=$(date +%s.%N)
awk "BEGIN {printf \"MovingSurface runtime [seconds]: %.5f\n\", $res2-$res1}"

# i defines how many free surface evolution time steps shall be performed
for i in {1..30}
do
    i2=$(printf %04d $i)

    # Generate updated 3D mesh
    res1=$(date +%s.%N)
    gmsh-aj -2 -format stl $name_dir$"surface_shell_domain.geo"
    python make_mesh.py
    echo $name_dir$"surface_2D_"$i2$".geo"
    gmsh-aj -2 $name_dir$"surface_2D_"$i2$".geo"
    gmsh-aj -3 $name_dir$"vol_mesh_domain.geo"
    dolfin-convert $name_dir$"vol_mesh_domain.msh" $name_dir$name$i2$name2
    dolfin-convert $name_dir$"surface_2D_"$i2$".msh" $name_dir$"surface_2D_"$i2$".xml"
    res2=$(date +%s.%N)
    awk "BEGIN {printf \"Mesh runtime [seconds]: %.5f\n\", $res2-$res1}"
    python read_mesh.py

    # Compute 3D velocity field
    res1=$(date +%s.%N)
    mpirun -np 4 python ice_dynamics3D.py
    rm $name_dir$"surface_shell_domain.stl"
    rm $name_dir$"surface_2D_"$i2$".geo"
    k=$(expr $i + 1)
    k2=$(printf %04d $k)
    cp surface_start1_shell.geo $name_dir$"surface_shell_domain.geo"
    cp surface_start1_2D.geo $name_dir$"surface_2D_"$k2$".geo"
    echo "cp new file"
    echo "cp surface_start1_2D.geo" $name_dir$"surface_2D_"$k2$".geo"
    res2=$(date +%s.%N)
    awk "BEGIN {printf \"Iceflow runtime [seconds]: %.5f\n\", $res2-$res1}"

    # Compute free surface evolution
    res1=$(date +%s.%N)
    python moving_surface.py
    res2=$(date +%s.%N)
    awk "BEGIN {printf \"MovingSurface runtime [seconds]: %.5f\n\", $res2-$res1}"

done

res4=$(date +%s.%N)
awk "BEGIN {printf \"Overall runtime [seconds]: %.5f\n\", $res4-$res3}"

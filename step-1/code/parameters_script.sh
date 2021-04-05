#!/bin/sh
# Refresh parameters

#change parameters in mesh_script.geo
Height=$(sed -n 1p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Height =/c Height = $Height;" data/mesh_script.geo
Radius=$(sed -n 2p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius =/c Radius = $Radius;" data/mesh_script.geo
Radius_Inner=$(sed -n 3p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius_Inner =/c Radius_Inner = $Radius_Inner;" data/mesh_script.geo
Nb=$(sed -n 4p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Nb =/c Nb = $Nb;" data/mesh_script.geo
Nc1=$(sed -n 5p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Nc1 =/c Nc1 = $Nc1;" data/mesh_script.geo
Nc2=$(sed -n 6p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Nc2 =/c Nc2 = $Nc2;" data/mesh_script.geo
Nz=$(sed -n 7p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Nz =/c Nz = $Nz;" data/mesh_script.geo

#Write Parameters in main
sed -i "/double height =/c double height = $Height;" code/main.cpp
sed -i "/double int_rad =/c double int_rad = $Radius_Inner;" code/main.cpp
sed -i "/double out_rad =/c double out_rad = $Radius;" code/main.cpp

echo -e 'Done \n'

#!/bin/sh
# Refresh parameters

#change parameters in mesh_script.geo
Height=$(sed -n 1p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Height =/c Height = $Height" data/mesh_script.geo

Radius=$(sed -n 2p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius =/c Radius = $Radius" data/mesh_script.geo

Radius_Inner=$(sed -n 3p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius_Inner =/c Radius_Inner = $Radius_Inner" data/mesh_script.geo

echo $Height
echo $Radius
echo $Radius_Inner
echo -e 'Done \n'

#!/bin/sh
# Refresh parameters

# Change parameters in mesh_script.geo
Height=$(sed -n 3p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Height =/c Height = $Height;" script/mesh_script.geo
Radius=$(sed -n 4p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius =/c Radius = $Radius;" script/mesh_script.geo
Radius_Inner=$(sed -n 5p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius_Inner =/c Radius_Inner = $Radius_Inner;" script/mesh_script.geo

Nb=$(sed -n 6p data/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nb =/c Nb = $Nb;" script/mesh_script.geo
Nc1=$(sed -n 7p data/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nc1 =/c Nc1 = $Nc1;" script/mesh_script.geo
Nc2=$(sed -n 7p data/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nc2 =/c Nc2 = $Nc2;" script/mesh_script.geo
Nz=$(sed -n 8p data/parameters.txt | tr -d -c 0-9.)
sed -i "/Nz =/c Nz = $Nz;" script/mesh_script.geo

echo -e 'Done! \n'

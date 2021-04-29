#!/bin/sh
# Refresh parameters

# Change parameters in mesh_script.geo
Height=$(sed -n 2p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/Height =/c Height = $Height;" settings/mesh_script.geo
Radius=$(sed -n 3p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius =/c Radius = $Radius;" settings/mesh_script.geo
Radius_Inner=$(sed -n 4p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/Radius_Inner =/c Radius_Inner = $Radius_Inner;" settings/mesh_script.geo

Nb=$(sed -n 5p settings/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nb =/c Nb = $Nb;" settings/mesh_script.geo
Nc1=$(sed -n 6p settings/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nc1 =/c Nc1 = $Nc1;" settings/mesh_script.geo
Nc2=$(sed -n 6p settings/parameters.txt | tr -d -c 0-9.)+1
sed -i "/Nc2 =/c Nc2 = $Nc2;" settings/mesh_script.geo
Nz=$(sed -n 7p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/Nz =/c Nz = $Nz;" settings/mesh_script.geo

echo -e 'Done! \n'

#/bin/sh
# Refresh parameters

# Change parameters in mesh_script.geo
Mesh_script=$(sed -n 2p settings/parameters.txt | cut -d ' ' -f 1)

R=$(sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Rmax =/c Rmax = $R;" $Mesh_script
Z=$(sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Zmax =/c Zmax = $Z;" $Mesh_script

l=$(sed -n 5p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/l =/c l = $l;" $Mesh_script
h=$(sed -n 6p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/h =/c h = $h;" $Mesh_script

NR=$(sed -n 7p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/NR =/c NR = $NR;" $Mesh_script
NZ=$(sed -n 8p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/NZ =/c NZ = $NZ;" $Mesh_script

echo -e 'Done! \n'

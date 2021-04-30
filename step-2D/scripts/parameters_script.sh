#/bin/sh
# Refresh parameters

# Change parameters in mesh_script.geo
Mesh_script=$(sed -n 2p data/parameters.txt | cut -d ' ' -f 1)
Rmin=$(sed -n 3p data/parameters.txt | tr -d -c 0-9.-)
sed -i "/Rmin =/c Rmin = $Rmin;" $Mesh_script
Rmax=$(sed -n 4p data/parameters.txt | tr -d -c 0-9.-)
sed -i "/Rmax =/c Rmax= $Rmax;" $Mesh_script
Zmin=$(sed -n 5p data/parameters.txt | tr -d -c 0-9.-)
sed -i "/Zmin =/c Zmin = $Zmin;" $Mesh_script
Zmax=$(sed -n 6p data/parameters.txt | tr -d -c 0-9.-)
sed -i "/Zmax =/c Zmax = $Zmax;" $Mesh_script

NR=$(sed -n 7p data/parameters.txt | tr -d -c 0-9.)
sed -i "/NR =/c NR = $NR;" $Mesh_script
NZ=$(sed -n 8p data/parameters.txt | tr -d -c 0-9.)
sed -i "/NZ =/c NZ = $NZ;" $Mesh_script

echo -e 'Done! \n'

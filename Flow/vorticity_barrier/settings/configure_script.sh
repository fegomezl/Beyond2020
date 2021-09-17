#/bin/sh
# Refresh parameters

# Change parameters in mesh_script.geo
Mesh_script=$(sed -n 2p settings/parameters.txt | cut -d ' ' -f 1)
Rmin=$(sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Rmin =/c Rmin = $Rmin;" $Mesh_script
Rmax=$(sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Rmax =/c Rmax = $Rmax;" $Mesh_script
Zmin=$(sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Zmin =/c Zmin = $Zmin;" $Mesh_script
Zmax=$(sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Zmax =/c Zmax = $Zmax;" $Mesh_script
Hside=$(sed -n 7p settings/parameters.txt | tr -d -c 0-9.-)
sed -i "/Hside =/c Hside = $Hside;" $Mesh_script

NR=$(sed -n 8p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/NR =/c NR = $NR;" $Mesh_script
NZ=$(sed -n 9p settings/parameters.txt | tr -d -c 0-9.)
sed -i "/NZ =/c NZ = $NZ;" $Mesh_script

echo -e 'Done! \n'

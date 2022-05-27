#/bin/sh

# Change parameters in mesh_script.geo
Mesh_script=$(sed -n 11p settings/parameters.txt | cut -d ' ' -f 1)
R=$(sed -n 12p settings/parameters.txt | tr -d -c 0-9.-)
Z=$(sed -n 13p settings/parameters.txt | tr -d -c 0-9.-)
r=$(sed -n 14p settings/parameters.txt | tr -d -c 0-9.)

R=$((R/r))
Z=$((Z/r))
r=$((r/r))

sed -i "/R =/c R = $R;" $Mesh_script
sed -i "/Z =/c Z = $Z;" $Mesh_script
sed -i "/r =/c r = $r;" $Mesh_script

echo -e 'Done! \n'

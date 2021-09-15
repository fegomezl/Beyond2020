#!/bin/bash

##Mesh parmeters
SCRIPT=$(sed -n 2p settings/parameters.txt | cut -d '#' -f 1)
RMIN=$(sed -n 3p settings/parameters.txt | tr -d -c 0-9.-)
RMAX=$(sed -n 4p settings/parameters.txt | tr -d -c 0-9.-)
ZMIN=$(sed -n 5p settings/parameters.txt | tr -d -c 0-9.-)
ZMAX=$(sed -n 6p settings/parameters.txt | tr -d -c 0-9.-)

##Program parameters
REF_I=$(sed -n 11p settings/parameters.txt | tr -d -c 0-9.)
REF_L=$(sed -n 12p settings/parameters.txt | tr -d -c 0-9.)
ORDER=$(sed -n 13p settings/parameters.txt | tr -d -c 0-9.)
T_FU=$(sed -n 16p settings/parameters.txt | tr -d -c 0-9.-)
C_L=$(sed -n 17p settings/parameters.txt | tr -d -c 0-9.)
C_S=$(sed -n 18p settings/parameters.txt | tr -d -c 0-9.)
K_L=$(sed -n 19p settings/parameters.txt | tr -d -c 0-9.)
K_S=$(sed -n 20p settings/parameters.txt | tr -d -c 0-9.)
L=$(sed -n 21p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_I=$(sed -n 22p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_JS=$(sed -n 23p settings/parameters.txt | tr -d -c 0-9.)
DELTA_T_J=$(sed -n 24p settings/parameters.txt | tr -d -c 0-9.)
DT=$(sed -n 25p settings/parameters.txt | tr -d -c 0-9.)
T_FI=$(sed -n 26p settings/parameters.txt | tr -d -c 0-9.)
VIS=$(sed -n 27p settings/parameters.txt | tr -d -c 0-9.)
ODE=$(sed -n 28p settings/parameters.txt | tr -d -c 0-9.)
RELT=$(sed -n 29p settings/parameters.txt | tr -d -c 0-9.)
ABST=$(sed -n 30p settings/parameters.txt | tr -d -c 0-9.)
STATE=$(sed -n 31p settings/parameters.txt | tr -d -c 0-9.)

##Compiling parameters
mpiexec -np ${PROCCESORS} ./main.x --mesh results/mesh.msh -Rmin ${RMIN} -Rmax ${RMAX} -Zmin ${ZMIN} -Zmax ${ZMAX} \
    -o ${ORDER} -r_i ${REF_I} -r_l ${REF_L} -DT_i ${DELTA_T_I} -DT_js ${DELTA_T_JS} -DT_j ${DELTA_T_J} \ 
    -T_f ${T_FU} -c_s ${C_S} -c_l ${C_L} -k_s ${K_S} -k_l ${K_L} -L ${L} \
    -dt ${DT} -t_f ${T_FI} -v_s ${VIS} -ode ${ODE} -reltol ${RELT} -abstol ${ABST} -state ${STATE}

##aca eestamos asumiendo wue el directorio de trabajo es ~/
##para trabajar en otro directorio hay que editar la direccion de los ejecutables, el mesh
## la bandera -cwd quiere decir, trabajar en el directorio donde esta el script
#!/bin/bash
#$ -pe mpi 16
#$ -q normal.q
#$ -cwd
#$ -S /bin/bash
#$ -o ./file.out
#$ -e ./file.err
##EJ.####
#mpiexec -v -np $NSLOTS  --oversubscribe -host hercules4,hercules3,hercules2 ./pi

mpiexec -np $NSLOTS --oversubscribe -host hercules4,hercules3,hercules2 \
~/brinicle/step-2/onephase_3/./main.x --mesh ~/brinicle/step-2/onephase_3/results/mesh.msh \
-Rmin 0 -Rmax 10 -Zmin 0 -Zmax 10 -o ${order} -r ${refinements} -T_f -10 -a_l 7.8 -a_s 70.8 \
-dt 0.01 -t_f 1 -v_s 10 -ode ${solver} -reltol 0.00001 -abstol 0.00001


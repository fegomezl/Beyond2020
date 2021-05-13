##aca eestamos asumiendo wue el directorio de trabajo es ~/  
##para trabajar en otro directorio hay que editar la direccion de los ejecutables, el mesh
## la bandera -cwd quiere decir, trabajar en el directorio donde esta el script
##la ejecucion desde la linea de comandos es -> qsub -v order=1,refinements=1,solver=12 run.sh
#!/bin/bash
#$ -pe mpi 10
#$ -q normal.q@hercules1
#$ -cwd
#$ -S /bin/bash
#$ -o ./file.out
#$ -e ./file.err
##EJ.####
##mpiexec -v -np $NSLOTS  --oversubscribe -host hercules4,hercules3,hercules2 ./pi

mpiexec -v -np $NSLOTS --oversubscribe -host hercules1\
	~/brinicle/step-1/cluster_3D/./main.x --mesh ~/brinicle/step-1/cluster_3D/results/mesh.msh --internal-radius 5 --outer-radius 10 -he 10 -o 1 -r 5

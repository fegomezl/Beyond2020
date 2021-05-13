#!/bin/bash
#$ -pe mpi 8
#$ -q normal.q
#$ -cwd
#$ -S /bin/bash
#$ -o ./file.out
#$ -e ./file.err
mpiexec -v -np $NSLOTS --oversubscribe -host hercules4,hercules3,hercules2 \
~/brinicle/step-1/cluster_3D/./main.x --mesh ~/brinicle/step-1/cluster_3D/results/mesh.msh \
--height 10   --internal-radius 5 --outer-radius 10 --order 2 --refinements 1

#!/bin/bash

for SOLVER in 8 9 11 12; do
		for REFINEMENT in 0 1 2 3 4 5; do
			qsub -v order=1,refinements=$REFINEMENT,solver=$SOLVER  run.sh
		done
done
#asi como esta todo va a quedar gurdado en un archivo llamado file.out, un run despues del otro
#para definir el nombre del archovo se usa -o, por ejemplo -o file.$ORDER.$REFINEMENT.$SOLVER
#para que esto funcione en run.sh hay que añadir un asterisco mas a la linea #$ -o file.out
for SOLVER in 8 9 11 12; do
        for ORDER in 1 2 3 4 5; do
                        qsub -v order=$ORDER,refinements=2,solver=$SOLVER  run.sh
                done
done

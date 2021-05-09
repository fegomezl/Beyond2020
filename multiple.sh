#!/bin/bash

for SOLVER in 12; do
	for ORDER in  2; do
		for REFINEMENT in 3 4; do
			qsub -v order=$ORDER,refinements=$REFINEMENT,solver=$SOLVER  run.sh
			done
		done
done
#asi como esta todo va a quedar gurdado en un archivo llamado file.out, un run despues del otro
#para definir el nombre del archovo se usa -o, por ejemplo -o file.$ORDER.$REFINEMENT.$SOLVER
#para que esto funcione en run.sh hay que añadir un asterisco mas a la linea #$ -o file.out

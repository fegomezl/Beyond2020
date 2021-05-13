#!/bin/bash

for ORDER in  1 2 3; do
	qsub -v order=$ORDER,refinements=4 run.sh
done
#asi como esta todo va a quedar gurdado en un archivo llamado file.out, un run despues del otro
#para definir el nombre del archovo se usa -o, por ejemplo -o file.$ORDER.$REFINEMENT.$SOLVER
#para que esto funcione en run.sh hay que añadir un asterisco mas a la linea #$ -o file.out

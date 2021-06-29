#!/bin/bash

PROCCESORS_O=16

qsub -cwd -S /bin/bash -q normal.q@hercules1 -N brinicle -o results/brinicle.out -e results/brinicle.err -pe make $PROCCESORS_O -v PROCCESORS=$PROCCESORS_O settings/job.sh &


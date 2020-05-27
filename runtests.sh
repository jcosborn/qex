#!/bin/sh

export OMP_NUM_THREADS=2
export RUNJOB="mpiexec -n 2 -bind-to user:0,1"
#export RUNJOB="mpiexec -n 2 --bind-to none"
#export RUNJOB="true"

./testscript.sh

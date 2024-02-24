#!/bin/bash
#$ -o ./log/$JOB_NAME_log.out

#$ -q openmp.q@phase08

#$ -pe mpi 1
#$ -j y
#$ -cwd

data=$1

. ~/.bashrc

mpirun -np 1 ./Lcoeff "${data}"

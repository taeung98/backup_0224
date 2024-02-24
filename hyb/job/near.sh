#!/bin/bash
#$ -o ./log/$JOB_NAME_log.out

#$ -q mpi.q@phase01

# we can use (openmp.q@phase05~09) or (mpi.q@phase01~07),NOTE that same phase0# is same cpu
# useful command
# qstat -u '*', qstf(more simple) : all user using q can be checked
# input qconf -sq (openmp.q)or(mpi.q) and cherck 'slots' that we are available

#$ -pe mpi 1 
#$ -j y
#$ -cwd

t0=$(date +%s.%N)
t0_string=$(date)

Nb=$1
near=$2
i1=$3
i2=$4
i3=$5

. ~/.bashrc
mpirun -np 1 ./hyb "$Nb" "$near" "$i1" "$i2" "$i3"

#########################################

t1=$(date +%s.%N)
t1_string=$(date)

t=$(echo "$t1 - $t0"|bc)
h=$(echo "($t/3600)"|bc)
m=$(echo "($t%3600)/60"|bc)
s=$(echo "($t%3600)%60"|bc)

echo ""
echo "# Job Start : $t0_string"
echo "# Job End   : $t1_string"
echo "# Elapsed time : ${h}h ${m}m ${s}s"
echo ""

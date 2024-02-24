#!/bin/bash

queue=$1-1
Nb=$2

for x in {1..3}; do
	for i in {1..6}; do
		for j in {1..6}; do
			for k in {1..2}; do

				qsub_cmd="qsub -q mpi.q@phase0$((${x}+${queue})) -N L_Nb${Nb}_n4_${i}${j}${k} job/Lcoeff.sh database_test/Nb${Nb}/Nb${Nb}_n{4}_${i}${j}${k}.h5"

				echo ${qsub_cmd}
#			${qsub_cmd}

			done
		done
	done
done

for x in {1..18}; do
	for y in {1..4};do

		qsub_cmd="qsub -q mpi.q@phase0$((${queue}+4)) -N L_Nb${Nb}_n123_${x}${y} job/Lcoeff.sh database_test/Nb${Nb}/n123_${x}${y}.h5"

#		echo ${qsub_cmd}
		${qsub_cmd}

	done
done

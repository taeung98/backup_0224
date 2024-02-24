#!/bin/bash

queue=$1-1
Nb=$2
near=$3

if [ $near == 4 ]; then

	for x in {1..3}; do
		for i in {1..6}; do
			for j in {1..6}; do
				for k in {1..2}; do

#					if [ $exe == "n" ]; then

					qsub_cmd="qsub -q mpi.q@phase0$((${x}+${queue})) -N Nb${Nb}_n${near}_${i}${j}$((2 * ${x} - (2 - ${k}))) job/near.sh ${Nb} ${near} ${i} ${j} $((2 * ${x} - (2 - ${k})))"

					echo ${qsub_cmd}
#					${qsub_cmd}

#					fi

				done
			done
		done
	done

else
	for x in {1..18}; do
		for y in {1..4};do
	
			qsub_cmd="qsub -q openmp.q@phase0$((${queue}+1)) -N Nb${Nb}_n${near}_${x}${y} job/near.sh ${Nb} ${near} ${x} ${y} 0"

#			echo ${qsub_cmd}
			${qsub_cmd}
		done
	done

fi

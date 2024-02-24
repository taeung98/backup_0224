#!/bin/bash

for i in {1..3}; do 
	for j in {1..2}; do 
		echo $((${i} * 2 - (2 - ${j})));
	done;
done

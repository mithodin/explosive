#!/bin/bash

for group in $(h5ls simulation_data.h5 | cut -f 1 -d " ")
do
	i=$(echo $group | rev | cut -d "-" -f 1 | rev)
	cd hy_$i
	h52ascii ../simulation_data.h5 /$group
	h5clustersize ../simulation_data.h5 /$group 10
	cd ..
	echo $group
done

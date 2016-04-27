#!/bin/bash
if [ $# -ne 3 ]
then
	echo "Usage: run_hysteresis <dirname> <template> <squarewell/continuous>"
else
	root=$(pwd)"/"
	for conffile in $(ls $2*)
	do
		sim_rundir="/home/lucas/Simulation/run/$1/$conffile/"
		mkdir -p $sim_rundir
		echo $sim_rundir
		cd /home/lucas/Simulation/bin_$3/build
		cp $root$conffile ../config.h
		make clean
		make -j
		cp ../bin/explosive $sim_rundir
		cd $sim_rundir
		./explosive 1> log 2> error.log
	done
fi

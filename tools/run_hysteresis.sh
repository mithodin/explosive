#!/bin/bash
root=$(pwd)"/"
for conffile in "$@"
do
	sim_rundir="/home/lucas/Simulation/run/hysteresis/$conffile/"
	mkdir -p $sim_rundir
	echo $sim_rundir
	cd /home/lucas/Simulation/engine/build
	cp $root$conffile ../config.h
	make -j
	cp ../bin/explosive $sim_rundir
	cd $sim_rundir
	./explosive 1> log 2> error.log
done

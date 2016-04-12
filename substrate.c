#include <x86intrin.h>
#include <stdbool.h>
#include "config.h"
#include "colloid.h"
#include "geometry.h"

vector2d centre;

int external_energy(vector2d position){
	double d=0;
	distance(position,centre,&d);
	if(d<10.0){ return -1; }
	return 0;
}

void init_substrate(void){
	centre[0]=SIZE_X/2.0;
	centre[1]=SIZE_Y/2.0;
}

/** @file */
#include <x86intrin.h>
#include <stdbool.h>
#include "config.h"
#include "colloid.h"
#include "geometry.h"

vector2d centre; /**< the centre of the well */

/**
 * Calculate the external energy of a colloid
 * 
 * @param position The 2d position of the colloid
 * @return The energy in units of well depth (ENERGY_WELL_DEPTH)
 */
int external_energy(vector2d position){
	double d=0;
	distance(position,centre,&d);
	if(d<10.0){ return -1; }
	return 0;
}

/**
 * Initialize the substrate
 */
void init_substrate(void){
	centre[0]=SIZE_X/2.0;
	centre[1]=SIZE_Y/2.0;
}

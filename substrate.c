/** @file */
#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "config.h"
#include "colloid.h"
#include "geometry.h"

vector2d well[SUBSTRATE_WELLS_X][SUBSTRATE_WELLS_Y]; /**< the centre of the wells */

/**
 * Calculate the external energy of a colloid
 * 
 * @param position The 2d position of the colloid
 * @return The energy in units of well depth (ENERGY_WELL_DEPTH)
 */
int external_energy(vector2d position){
	int row,column;
	row=((int)(position[1]*SUBSTRATE_WELLS_Y/SIZE_Y))%SUBSTRATE_WELLS_Y;
	column=((int)(position[0]*SUBSTRATE_WELLS_X/SIZE_X+(1.0-(row%2))/2.0))%SUBSTRATE_WELLS_X;
	double d=0;
	distance(position,well[column][row],&d);
	if(d<SUBSTRATE_WELL_RADIUS){ return -1; }
	return 0;
}

/**
 * Initialize the substrate
 */
void init_substrate(void){
	for(int row=0;row<SUBSTRATE_WELLS_Y;++row){
		for(int column=0;column<SUBSTRATE_WELLS_X;++column){
			well[column][row][0]=((row%2)*SUBSTRATE_OFFSET_ODD)+column*(SIZE_X/SUBSTRATE_WELLS_X);
			well[column][row][1]=(row+0.5)*(SIZE_Y/SUBSTRATE_WELLS_Y);
		}
	}
}

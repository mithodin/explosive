/** @file */
#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include <stdio.h>
#include "config.h"
#include "geometry.h"
#include "colloid.h"

vector2d well[SUBSTRATE_WELLS_X][SUBSTRATE_WELLS_Y]; /**< the centre of the wells */
vector4d unit_vectors[2];
double normalize;

void sample_energy(void);

/**
 * Calculate the external energy of a colloid
 * 
 * @param position The 2d position of the colloid
 * @return The energy in units of well depth (ENERGY_WELL_DEPTH)
 */
double external_energy(vector2d position){
	vector4d position2=_mm256_broadcast_pd(&(position.v));
	vector4d p1=_mm256_mul_pd(position2,unit_vectors[0]),p2=_mm256_mul_pd(position2,unit_vectors[1]);
	vector4d sum=_mm256_hadd_pd(p1,p2);
	vector4d pcos=_mm256_cos_pd(sum);
	vector4d hsum=_mm256_hadd_pd(pcos,pcos);
	double energy=(hsum[0]+hsum[2])/(-4.0);
	if( energy < 0.0 ){ energy*= normalize; }
	return energy<-ENERGY_WELL_DEPTH?-ENERGY_WELL_DEPTH:energy;
}

/**
 * Initialize the substrate
 */
void init_substrate(void){
	normalize=1.0;
	double scale=4.0*M_PI*SUBSTRATE_WELLS_X/SIZE_X/sqrt(3.0);
	double u0[4]={0.0,scale,sqrt(3.0)/2.0*scale,scale/2.0};
	double u1[4]={sqrt(3.0)/2.0*scale,-scale/2.0,0.0,0.0};
	unit_vectors[0]=_mm256_load_pd(u0);
	unit_vectors[1]=_mm256_load_pd(u1);

	vector2d test={_mm_set_pd(0.0,SUBSTRATE_WELL_RADIUS)};
	normalize=-ENERGY_WELL_DEPTH/external_energy(test);
	sample_energy();
}

void sample_energy(void){
	FILE *energy_file = fopen("potential.dat","w");
	vector2d test;
	for(int i=0;i<100;++i){
		test.c.x = i*SIZE_X/100.0/SUBSTRATE_WELLS_X;
		for(int j=0;j<200;++j){
			test.c.y = j*SIZE_Y/100.0/SUBSTRATE_WELLS_Y;
			fprintf(energy_file,"%1.5f\t%1.5f\t%1.10f\n",test.c.x,test.c.y,external_energy(test));
		}
	}
}

/*
set xrange [0:100]
set yrange [0:86.60254037844386467635]
size_x=100.0
num_x=5
scale=2.0/sqrt(3.0)*2.0*pi*num_x/size_x
g1x=0.0*scale
g1y=1.0*scale
g2x=sqrt(3.0)/2.0*scale
g2y=0.5*scale
g3x=g2x
g3y=-g2y

g(x,y)=(cos(g1x*x+g1y*y)+cos(g2x*x+g2y*y)+cos(g3x*x+g3y*y)+1)/(-4.0)
set contour
unset surface
set isosamples 100
set view equal xy
set view 0,0
set xlabel "x"
set ylabel "y"

splot g(x,y)
*/

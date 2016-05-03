/** @file */
#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include <stdio.h>
#include "config.h"
#include "geometry.h"
#include "colloid.h"

#ifdef SUBSTRATE_TRIGONAL

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

#endif

#ifdef SUBSTRATE_RANDOM
#include <stdlib.h>
#include "dSFMT/dSFMT.h"
#include "globals.h"

bool substrate_collision(vector2d *, vector2d, int);
double energy_single_well(double);
double energy_raw(double);
void sample_energy(void);
double grid_interpolate(vector2d, vector4d, vector4d, vector4d);
double energy_substrate_direct(vector2d);

vector4d *coefficients, *points_x, *points_y;
int grid_res_x,grid_res_y;
double cutoff_r,dx,dy;
vector2d patches[SUBSTRATE_NUMBER_OF_PATCHES];

/**
 * Calculate the external energy of a colloid
 * 
 * @param position The 2d position of the colloid
 * @return The energy in units of well depth (ENERGY_WELL_DEPTH)
 */
double external_energy(vector2d position){
	int row, column, i;
	column = (int)(position.c.x/dx);
	row = (int)(position.c.y/dy);
	i = row*grid_res_x+column;
	return grid_interpolate(position, points_x[i], points_y[i], coefficients[i]);
}

/**
 * Initialize the substrate
 */
void init_substrate(void){
	printf("> initializing substrate randomly with %d patches\n", SUBSTRATE_NUMBER_OF_PATCHES);
	grid_res_x=(int)(ceil(SIZE_X/0.1));
	grid_res_y=(int)(ceil(SIZE_Y/0.1));

	cutoff_r=SIZE_X<SIZE_Y?SIZE_X/2.0:SIZE_Y/2.0;

	for(int i=0;i<SUBSTRATE_NUMBER_OF_PATCHES;++i){
		do{
			patches[i].v=_mm_set_pd(SIZE_Y*dsfmt_genrand_open_close(&rng),SIZE_X*dsfmt_genrand_open_close(&rng));
		}while(substrate_collision(patches, patches[i], i));
	}
	coefficients=calloc(grid_res_x*grid_res_y,sizeof(vector4d));
	points_x=calloc(grid_res_x*grid_res_y,sizeof(vector4d));
	points_y=calloc(grid_res_x*grid_res_y,sizeof(vector4d));
	vector2d r11,r12,r21,r22;
	dx=SIZE_X/grid_res_x;
	dy=SIZE_Y/grid_res_y;
	vector2d rx, ry;
	for(int i=0;i<grid_res_x*grid_res_y;++i){
		r11.c.x=(i%grid_res_x)*dx;
		r11.c.y=(i/grid_res_x)*dy;
		r12.v = _mm_set_pd(r11.c.y+dy,r11.c.x);
		r21.v = _mm_set_pd(r11.c.y,r11.c.x+dx);
		r22.v = _mm_set_pd(r11.c.y+dx,r11.c.x+dx);
		points_x[i]=_mm256_set_pd(r11.c.x,r21.c.x,r11.c.x,r21.c.x);
		points_y[i]=_mm256_set_pd(r11.c.y,r11.c.y,r12.c.y,r12.c.y);
		coefficients[i]=_mm256_set_pd(energy_substrate_direct(r22)/dx/dy,-energy_substrate_direct(r12)/dx/dy,-energy_substrate_direct(r21)/dx/dy,energy_substrate_direct(r11)/dx/dy);
	}
	sample_energy();
}

double grid_interpolate(vector2d r, vector4d x, vector4d y, vector4d coeff){
	vector4d xx=_mm256_broadcast_sd(&(r.c.x));
	vector4d yy=_mm256_broadcast_sd(&(r.c.y));
	xx=_mm256_sub_pd(x,xx);
	yy=_mm256_sub_pd(y,yy);
	xx=_mm256_mul_pd(xx,yy);
	xx=_mm256_mul_pd(xx,coeff);
	xx=_mm256_hadd_pd(xx,xx);
	return xx[0]+xx[2];
}

bool substrate_collision(vector2d *patches, vector2d new, int i){
	double d=0;
	for(i--;i>=0;--i){
		distance(new,patches[i],&d);
		if( d < 2.0*SUBSTRATE_WELL_RADIUS ){
			return true;
		}
	}
	return false;
}

double energy_substrate_direct(vector2d r){
	double d=0.0,e=0.0;
	for(int i=0;i<SUBSTRATE_NUMBER_OF_PATCHES;++i){
		distance(r,patches[i],&d);
		e+=energy_single_well(d);
	}
	return e;
}

double energy_single_well(double distance){
	double frc=energy_raw(cutoff_r);
	return distance<SUBSTRATE_WELL_RADIUS?-ENERGY_WELL_DEPTH:(distance>cutoff_r?0.0:(frc-energy_raw(distance))/(frc/ENERGY_WELL_DEPTH-1.0));
}

double energy_raw(double distance){
	return -0.04*ENERGY_WELL_DEPTH/(distance+0.2-SUBSTRATE_WELL_RADIUS)/(distance+0.2-SUBSTRATE_WELL_RADIUS);
}

void sample_energy(void){
	FILE *energy_file = fopen("potential.dat","w");
	vector2d test;
	for(int i=0;i<1000;++i){
		test.c.x = i*SIZE_X/1000.0;
		for(int j=0;j<1000;++j){
			test.c.y = j*SIZE_Y/1000.0;
			//fprintf(energy_file,"%1.5f\t%1.5f\t%1.10f\n",test.c.x,test.c.y,external_energy(test));
			fprintf(energy_file,"%1.5f\t%1.5f\t%1.10f\t%1.10f\n",test.c.x,test.c.y,energy_substrate_direct(test),external_energy(test));
		}
	}
	fclose(energy_file);
}
#endif

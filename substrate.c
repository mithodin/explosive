/** @file */
#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "dSFMT/dSFMT.h"
#include "geometry.h"
#include "colloid.h"
#include "globals.h"
#include "logger.h"

#ifdef CONTINUE
#include <hdf5.h>
#include <hdf5_hl.h>
extern hid_t conf_group;
#endif

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
#if SUBSTRATE_NUMBER_OF_PATCHES > 0
	int row, column, i;
	column = (int)(position.c.x/dx);
	row = (int)(position.c.y/dy);
	i = row*grid_res_x+column;
	return grid_interpolate(position, points_x[i], points_y[i], coefficients[i]);
#else
	return 0.0;
#endif
}

/**
 * Initialize the substrate
 */
bool init_substrate(void){
#if SUBSTRATE_NUMBER_OF_PATCHES > 0
	#if SUBSTRATE_PATTERN == 0
	printf("> initializing substrate randomly with %d patches\n", SUBSTRATE_NUMBER_OF_PATCHES);
	#elif SUBSTRATE_PATTERN == 1
	printf("> initializing substrate trigonal pattern with %d patches\n", SUBSTRATE_NUMBER_OF_PATCHES);
	#elif SUBSTRATE_PATTERN == 2
	printf("> initializing substrate square pattern with %d patches\n", SUBSTRATE_NUMBER_OF_PATCHES);
	#endif
	grid_res_x=(int)(ceil(SIZE_X/0.1));
	grid_res_y=(int)(ceil(SIZE_Y/0.1));

	cutoff_r=SIZE_X<SIZE_Y?SIZE_X/2.0:SIZE_Y/2.0;

#ifdef CONTINUE
	double buffer[SUBSTRATE_NUMBER_OF_PATCHES*2];
	herr_t status = H5LTread_dataset_double(conf_group,"random_patches",buffer);
	if( status < 0 ){
		printf("> Error reading patch locations\n");
		return false;
	}
	for(int i=0;i<SUBSTRATE_NUMBER_OF_PATCHES;++i){
		patches[i].v=_mm_load_pd(&(buffer[2*i]));
	}
#else
	#if SUBSTRATE_PATTERN > 0
	int patches_x = (int)sqrt(SUBSTRATE_NUMBER_OF_PATCHES), patches_y = SUBSTRATE_NUMBER_OF_PATCHES/patches_x;
	double sx = SIZE_X/patches_x, sy = SIZE_Y/patches_y;
	#endif
	#if SUBSTRATE_PATTERN == 0
	int j=0;
	#endif
	for(int i=0;i<SUBSTRATE_NUMBER_OF_PATCHES;++i){
	#if SUBSTRATE_PATTERN == 0
		do{
			patches[i].v = _mm_set_pd(SIZE_Y*dsfmt_genrand_open_close(&rng),SIZE_X*dsfmt_genrand_open_close(&rng));
			if(++j > 1000000){
				i=0;
				break;
			}
		}while(substrate_collision(patches, patches[i], i));
		j=0;
	#elif SUBSTRATE_PATTERN == 1
		patches[i].v = _mm_set_pd((i/patches_x+0.5)*sy,(i%patches_x+0.5*((i/patches_x)%2))*sx);
	#elif SUBSTRATE_PATTERN == 2
		patches[i].v = _mm_set_pd((i/patches_x+0.5)*sy,(i%patches_x+0.5)*sx);
	#endif
	}

#endif
	if( !log_substrate(patches) ){
		printf("> Error storing patch locations\n");
		return false;
	};
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
#endif //END IF NO PATCHES
	return true;
}

double grid_interpolate(vector2d r, vector4d x, vector4d y, vector4d coeff){
	vector4d xx=_mm256_broadcast_sd(&(r.c.x));
	vector4d yy=_mm256_broadcast_sd(&(r.c.y));
	xx=_mm256_sub_pd(x,xx);
	yy=_mm256_sub_pd(y,yy);
	xx=_mm256_mul_pd(xx,yy);
	xx=_mm256_mul_pd(xx,coeff);
	double *res;
	posix_memalign(&res,32*sizeof(void *),sizeof(double));
	_mm256_store_pd(res,_mm256_hadd_pd(xx,xx));
	return res[0]+res[2];
}

bool substrate_collision(vector2d *patches, vector2d new, int i){
	double d=0;
	for(i--;i>=0;--i){
		distance(new,patches[i],&d);
		if( d < 2.0*SUBSTRATE_WELL_RADIUS+2.0 ){
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
	#if SUBSTRATE_CONTINUOUS == 1
	double frc=energy_raw(cutoff_r);
	return distance<SUBSTRATE_WELL_RADIUS?-ENERGY_WELL_DEPTH:(distance>cutoff_r?0.0:(frc-energy_raw(distance))/(frc/ENERGY_WELL_DEPTH-1.0));
	#elif SUBSTRATE_CONTINUOUS == 0
	return distance < SUBSTRATE_WELL_RADIUS?-ENERGY_WELL_DEPTH:0.0;
	#endif
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

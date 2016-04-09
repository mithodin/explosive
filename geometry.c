/** @file */

#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "colloid.h"
#include "geometry.h"

/**
 * Calculate the distance between two points
 *
 * @return The distance between two points
 */
double distance_direct(double x1,double y1,double x2,double y2){
	double dx=x1-x2;
	#ifdef PERIODIC_X
	dx=dx>(SIZE_X/2)?dx-SIZE_X:dx;
	#endif
	double dy=y1-y2;
	#ifdef PERIODIC_Y
	dy=dy>(SIZE_Y/2)?dy-SIZE_Y:dy;
	#endif
	return sqrt(dx*dx+dy*dy);
}

/**
 * Calculate the distance between two points using SIMD instructions
 *
 * @param r1 Starting vector
 * @param r2 End vector
 * @param distance The length of the relative vector
 * @return The vector pointing from r1 to r2
 */
vector2d distance(vector2d r1, vector2d r2, double *distance){
	vector2d dist=_mm_sub_pd(r2, r1); //return vector pointing from r1 to r2
	#ifdef PERIODIC_X
	if(dist[0]>(SIZE_X/2)){
		dist[0]-=SIZE_X;
	}
	#endif
	#ifdef PERIODIC_Y
	if(dist[1]>(SIZE_Y/2)){
		dist[1]-=SIZE_Y;
	}
	#endif
	vector2d dist_squared=_mm_mul_pd(dist,dist);
	*distance=sqrt(dist_squared[0]+dist_squared[1]);
	return dist;
}

void make_periodic(vector2d *r){
	#ifdef PERIODIC_X
	if((*r)[0] < 0){ (*r)[0]=fmod((*r)[0],SIZE_X)+SIZE_X; }
	if((*r)[0] > SIZE_X){ (*r)[0]=fmod((*r)[0],SIZE_X); }
	#endif
	#ifdef PERIODIC_Y
	if((*r)[1] < 0){ (*r)[1]=fmod((*r)[1],SIZE_Y)+SIZE_Y; }
	if((*r)[1] > SIZE_Y){ (*r)[1]=fmod((*r)[1],SIZE_Y); }
	#endif
}

double angle_twopi(double angle){
	if (angle>0){
		return fmod(angle, 2.0*M_PI);
	}
	else {
		return fmod(angle, 2.0*M_PI)+2.0*M_PI;
	}
}

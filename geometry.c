/** @file */

#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "geometry.h"
#include "colloid.h"

/**
 * Calculate the distance between two points
 *
 * @return The distance between two points
 */
double distance_direct(double x1,double y1,double x2,double y2){
	double diff_x=x1-x2;
	#ifdef PERIODIC_X
	diff_x=diff_x>(SIZE_X/2)?diff_x-SIZE_X:diff_x;
	#endif
	double diff_y=y1-y2;
	#ifdef PERIODIC_Y
	diff_y=diff_y>(SIZE_Y/2)?diff_y-SIZE_Y:diff_y;
	#endif
	return sqrt(diff_x*diff_x+diff_y*diff_y);
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
	vector2d dist={_mm_sub_pd(r2.v, r1.v)}; //return vector pointing from r1 to r2
	#ifdef PERIODIC_X
	if(dist.c.x>(SIZE_X/2)){
		dist.c.x-=SIZE_X;
	}else if(dist.c.x<(SIZE_X/(-2))){
		dist.c.x+=SIZE_X;
	}
	#endif
	#ifdef PERIODIC_Y
	if(dist.c.y>(SIZE_Y/2)){
		dist.c.y-=SIZE_Y;
	}else if(dist.c.y<(SIZE_Y/(-2))){
		dist.c.y+=SIZE_Y;
	}
	#endif
	vector2d dist_squared={_mm_mul_pd(dist.v,dist.v)};
	*distance=sqrt(dist_squared.c.x+dist_squared.c.y);
	return dist;
}

/**
 * Constrain a vector to a box with periodic boundary conditions.
 *
 * @param r a 2d vector
 */
vector2d make_periodic(vector2d r){
	#ifdef PERIODIC_Y
	if(r.c.y < 0){ r.c.y=fmod(r.c.y,SIZE_Y)+SIZE_Y; }
	if(r.c.y > SIZE_Y){ r.c.y=fmod(r.c.y,SIZE_Y); }
	#endif
	#ifdef PERIODIC_X
	if(r.c.x < 0){ r.c.x=fmod(r.c.x,SIZE_X)+SIZE_X; }
	if(r.c.x > SIZE_X){ r.c.x=fmod(r.c.x,SIZE_X); }
	#endif
	return r;
}

/**
 * Make an angle fall into the interval [0;2pi]
 *
 * @param angle an angle in any interval
 * @return an angle in interval [0;2pi]
 */
double angle_twopi(double angle){
	if (angle>0){
		return fmod(angle, 2.0*M_PI);
	}
	else {
		return fmod(angle, 2.0*M_PI)+2.0*M_PI;
	}
}

/**
 * Calculate the y distance between two points. Handles periodicity.
 *
 * @param y0 start point
 * @param y1 end point
 * @return the (periodic) distance y1-y0
 */
double distance_y(double y0, double y1){
	double dy=y1-y0;
	#ifdef PERIODIC_Y
	if(dy>(SIZE_Y/2)){
		dy-=SIZE_Y;
	}else if(dy<(SIZE_Y/(-2))){
		dy+=SIZE_Y;
	}
	#endif
	return dy;
}

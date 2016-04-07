/** @file */

#include <math.h>
#include "distance.h"

/**
 * Calculate the distance between two points
 * 
 * @return The distance between two points
 */
double distance(double x1,double y1,double x2,double y2){
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

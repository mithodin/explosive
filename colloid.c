#include <math.h>
#include <stdbool.h>
#include <x86intrin.h>
#include "config.h"
#include "dSFMT/dSFMT.h"
#include "colloid.h"
#include "globals.h"
#include "geometry.h"

vector2d colloid_patch_site(Colloid *, int);
int colloid_closest_site(double, double);

bool colloid_bonded(Colloid *c1, Colloid *c2){
	double distance_abs=0;
	vector2d relative_vector=distance(c1->position,c2->position,&distance_abs);
	if(distance_abs > COLLOID_MIN_BONDING_DISTANCE || distance_abs < COLLOID_DIAMETER){
		return false;
	}
	else{
		double beta1=atan2(relative_vector[1],relative_vector[0]);
		double beta2=beta1+M_PI;
		int c1i=colloid_closest_site(c1->phi,beta1);
		int c2i=colloid_closest_site(c2->phi,beta2);
		distance(colloid_patch_site(c1,c1i),colloid_patch_site(c2,c2i),&distance_abs);
		if(distance_abs<COLLOID_PATCH_DIAMETER){
			return true;
		}else{
			return false;
		}
	}
}

int colloid_closest_site(double alpha, double beta){
	return (int)(angle_twopi(beta-alpha+ONE_THIRD_PI)/TWO_THIRDS_PI);
}

vector2d colloid_patch_site(Colloid *c, int site){
	vector2d position;
	position[0]=c->position[0]+COLLOID_DIAMETER/2.0*cos(c->phi+TWO_THIRDS_PI*site);
	position[1]=c->position[1]+COLLOID_DIAMETER/2.0*sin(c->phi+TWO_THIRDS_PI*site);
	return position;
}

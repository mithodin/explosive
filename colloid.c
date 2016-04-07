#include <math.h>
	#include <stdio.h>
#include <stdbool.h>
#include <x86intrin.h>
#include "dSFMT/dSFMT.h"
#include "colloid.h"
#include "parameters.h"
#include "distance.h"

vector2d colloid_patch_site(Colloid *, int);
int colloid_closest_site(double, double);

bool colloid_bonded(Colloid *c1, Colloid *c2){
	double distance_abs=0;
	vector2d relative_vector=distance(c1->position,c2->position,&distance_abs);
	if(distance_abs > colloid_min_bond_dist || distance_abs < colloid_diameter){
printf("%d\t%d\t",-1,-1);
		return false;
	}
	else{
		double beta1=atan2(relative_vector[1],relative_vector[0]);
		double beta2=beta1+M_PI;
		int c1i=colloid_closest_site(c1->phi,beta1);
		int c2i=colloid_closest_site(c2->phi,beta2);
printf("%d\t%d\t",c1i,c2i);
		distance(colloid_patch_site(c1,c1i),colloid_patch_site(c2,c2i),&distance_abs);
		if(distance_abs<colloid_patch_diameter){
			return true;
		}else{
			return false;
		}
	}
}

int colloid_closest_site(double alpha, double beta){
	return (int)(angle_twopi(beta-alpha+onethirdpi)/twothirdspi);
}

double angle_twopi(double angle){
	if (angle>0){
		return fmod(angle, 2.0*M_PI);
	}
	else {
		return fmod(angle, 2.0*M_PI)+2.0*M_PI;
	}
}

vector2d colloid_patch_site(Colloid *c, int site){
	vector2d position;
	position[0]=c->position[0]+colloid_diameter/2.0*cos(c->phi+twothirdspi*site);
	position[1]=c->position[1]+colloid_diameter/2.0*sin(c->phi+twothirdspi*site);
	return position;
}

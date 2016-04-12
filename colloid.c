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

/**
 * Calculate whether two colloids are bonded. Also checks for collisions between colloids.
 *
 * @param c1 The first colloid
 * @param c2 The second colloid
 * @param collision A pointer to a bool to signal a collision
 * @param site1 A pointer to an int to return c1's bonding site
 * @param site2 A pointer to an int to return c2's bonding site
 * @return a bool signaling whether or not the colloids are bonded
 */
bool colloid_bonded(Colloid *c1, Colloid *c2, bool *collision, int *site1, int *site2){
	double distance_abs=0;
	*collision=false;
	vector2d relative_vector=distance(c1->position,c2->position,&distance_abs);
	if(distance_abs > COLLOID_MIN_BONDING_DISTANCE){
		return false;
	}else if(distance_abs < COLLOID_DIAMETER){
		*collision=true;
		return false;
	}else{
		double beta1=atan2(relative_vector[1],relative_vector[0]);
		double beta2=beta1+M_PI;
		int c1i=colloid_closest_site(c1->phi,beta1);
		int c2i=colloid_closest_site(c2->phi,beta2);
		distance(colloid_patch_site(c1,c1i),colloid_patch_site(c2,c2i),&distance_abs);
		if(distance_abs<COLLOID_PATCH_DIAMETER){
			if(site1!=NULL){
				*site1=c1i;
			}
			if(site2!=NULL){
				*site2=c2i;
			}
			return true;
		}else{
			return false;
		}
	}
}

/**
 * Calculate the index of the site that is closest to another colloid.
 *
 * @param alpha The current rotational position of the colloid
 * @param beta The angle between the vector pointing from the current colloid to the bonding colloid and the x axis.
 * @return index of the patch closest to the bonding colloid
 */
int colloid_closest_site(double alpha, double beta){
	return (int)(angle_twopi(beta-alpha+ONE_THIRD_PI)/TWO_THIRDS_PI);
}

/**
 * Calculate the (x,y) position of a given bonding site on a colloid
 * 
 * @param c The parent colloid
 * @param site The index of the patch
 * @return A 2d vector giving the position of the patch
 */
vector2d colloid_patch_site(Colloid *c, int site){
	vector2d position;
	position[0]=c->position[0]+COLLOID_DIAMETER/2.0*cos(c->phi+TWO_THIRDS_PI*site);
	position[1]=c->position[1]+COLLOID_DIAMETER/2.0*sin(c->phi+TWO_THIRDS_PI*site);
	return position;
}

/**
 * Initialize the double-linked-list of all colloids, sorted by y coordinate
 */
void init_ysorted_list(void){
	Colloid *root=&(particles[0]);
	particles[0].above=NULL;
	particles[0].below=NULL;
	#ifdef PERIODIC_Y
	Colloid *top=root;
	Colloid *bottom=root;
	#endif
	for(int i=1;i<NUMBER_OF_PARTICLES;++i){
		insert_sorted_y(&(particles[i]),root);
		#ifdef PERIODIC_Y
		if(particles[i].above==NULL){
			top=&(particles[i]);
		}
		if(particles[i].below==NULL){
			bottom=&(particles[i]);
		}
		#endif
	}
	#ifdef PERIODIC_Y
	top->above=bottom;
	bottom->below=top;
	#endif
}

/**
 * Insert a new colloid into the sorted double-linked list.
 *
 * @param c The new colloid to insert
 * @param list The entry point into the list
 */
void insert_sorted_y(Colloid *c, Colloid *list){
	Colloid *current=list;
	Colloid *tmp;
	c->above=NULL;
	c->below=NULL;
	while(current->above != NULL && current->above->position[1] >= current->position[1] && current->position[1]<c->position[1]){
		current=current->above;
	}
	while(current->below != NULL && current->below->position[1] <= current->position[1] && current->position[1]>c->position[1]){
		current=current->below;
	}
	if(current->position[1]>c->position[1]){
		tmp=current->below;
		current->below=c;
		c->above=current;
		c->below=tmp;
		if(tmp != NULL){
			tmp->above=c;
		}
	}else{
		tmp=current->above;
		current->above=c;
		c->below=current;
		c->above=tmp;
		if(tmp != NULL){
			tmp->below=c;
		}
	}
}

/**
 * Initialize all bonding sites on all colloids
 */
void init_bonding_partners(void){
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		particles[i].internal_energy=0;
		for(int k=0;k<3;++k){
			particles[i].bonding_partner[k]=NULL;
			particles[i].bond_site[k]=-1;
		}
	}
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		for(int j=i+1;j<NUMBER_OF_PARTICLES;++j){
			int site1=0;
			int site2=0;
			bool collision=false;
			if(colloid_bonded(&(particles[i]),&(particles[j]),&collision,&site1,&site2) && !collision){
				particles[i].bonding_partner[site1]=&(particles[j]);
				particles[j].bonding_partner[site2]=&(particles[i]);
				particles[i].bond_site[site1]=site2;
				particles[j].bond_site[site2]=site1;
				particles[i].internal_energy--;
				particles[j].internal_energy--;
			}
		}
	}
}

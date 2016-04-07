#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "colloid.h"
#include "monte_carlo.h"
#include "dSFMT/dSFMT.h"

int get_random_seed(void);

const double colloid_diameter=1.0;
const double colloid_patch_diameter=0.11965683746373795115;
const double colloid_min_bond_dist=1.11965683746373795115;
const double energy_well_depth=0.5;
const double energy_bond=1.0;
const double onethirdpi=1.04719755119659763132;
const double twothirdspi=2.09439510239319526264;
const int num_particles=2;

extern Colloid *particles;

dsfmt_t rng; /**< stores the RNG state */

int main(void){
	dsfmt_init_gen_rand(&rng,get_random_seed()); //initialize the rng

	mc_init(1.0);
	bool bonded=false;
	while(!bonded){
		mc_run(1);
		bonded=colloid_bonded(&particles[0],&particles[1]);
		printf("%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%1.5f\t%d\n",particles[0].position[0],particles[0].position[1],particles[0].phi,particles[1].position[0],particles[1].position[1],particles[1].phi,bonded?1:0);
	}
	return 0;
}

/**
 ** Get a random int to use as seed
 **
 ** If /dev/random is not available or fails to read, uses timeofday to get a "random" number
 **
 ** @return Returns a more or less random int
 **/
int get_random_seed(void){
	FILE *rnd = fopen("/dev/random","r");
	int seed;
	if( rnd ){
		if(fread(&seed,sizeof(int),1,rnd) != 1){
			struct timeval t;
			gettimeofday(&t,NULL);
			seed=((int)t.tv_usec)^((int)t.tv_sec);
		}
	}else{
		struct timeval t;
		gettimeofday(&t,NULL);
		seed=((int)t.tv_usec)^((int)t.tv_sec);
	}
	return seed;
}


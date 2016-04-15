/**
 * @file explosive.c
 * @brief The main program
 */

#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "config.h"
#include "colloid.h"
#include "logger.h"
#include "monte_carlo.h"
#include "dSFMT/dSFMT.h"

int get_random_seed(void);

//globals
Colloid particles[NUMBER_OF_PARTICLES]; /**< stores all the particles */
dsfmt_t rng; /**< stores the RNG state */
struct timeval sim_start_time; /**< stores the start of the simulation */

int main(void){
	gettimeofday(&sim_start_time,NULL);
	dsfmt_init_gen_rand(&rng,get_random_seed()); //initialize the rng

	double kbt=TEMPERATURE;
	double paccept;
	if( !mc_init(kbt) ){ printf("> mc_init() failed\n"); return -1; }
	if( !log_init() ){ printf("> log_init() failed.\n"); return -1; }
	paccept=mc_run(MONTE_CARLO_STEPS_MAIN,true);

	struct timeval sim_end_time;
	gettimeofday(&sim_end_time,NULL);
	log_simulation_stats((unsigned long)(sim_end_time.tv_sec-sim_start_time.tv_sec), paccept);
	if( !log_close()){ printf("> log_close() failed.\n"); return -1; }
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


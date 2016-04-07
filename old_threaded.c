/** @file */

#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "config.h"
#include "colloid.h"
#include "statistics.h"
#include "parameters.h"
#include "monte_carlo.h"
#include "initialize.h"
#import "dSFMT/dSFMT.h"

int getRandomSeed(void);

bool initdmax = true; /**< initialize amax and dmax? */

/**
 * Spawns a new thread with given parameters
 *
 * @param params Pointer to Config struct. Make sure to initialize correctly (see load_config.c)
 */
void *newThread(void *params){
	dsfmt_init_gen_rand(&(c->myrand),getRandomSeed()); //initialize the rng
}

/**
 * Get a random int to use as seed
 *
 * If /dev/random is not available or fails to read, uses timeofday to get a "random" number
 *
 * @return Returns a more or less random int
 */
int getRandomSeed(void){
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

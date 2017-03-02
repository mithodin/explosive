/**
 * @file globals.h
 * @brief All global variables are defined here
 */

extern dsfmt_t rng;
extern Colloid particles[NUMBER_OF_PARTICLES];
extern struct timeval sim_start_time;
extern int cluster_sizes[NUMBER_OF_PARTICLES];

extern int debug_int_energy;

//use these so we don't calculate them millions of times over
#define ONE_THIRD_PI 1.04719755119659763132
#define TWO_THIRDS_PI 2.09439510239319526264

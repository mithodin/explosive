#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "colloid.h"
#include "logger.h"
#include "globals.h"
#include "monte_carlo.h"
#include "geometry.h"

double monte_carlo_step(double, double);
bool mc_energy_change(vector2d, double, int, int*, int*);
double mc_acceptance_probability(int, int);
void mc_init_particles(void);
void mc_init_acceptance_probabilities(double);

double acceptance_probabilities_bonds[7];
double acceptance_probabilities_well[3];

double monte_carlo_step(double max_displacement, double max_rotation){
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		vector2d new_position;
		new_position[0]=particles[i].position[0]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		#ifndef PERIODIC_X
		if(new_position[0] < 0 || new_position[0] > SIZE_X){
			continue;
		}
		#endif
		new_position[1]=particles[i].position[1]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		#ifndef PERIODIC_Y
		if(new_position[1] < 0 || new_position[1] > SIZE_Y){
			continue;
		}
		#endif
		#if defined(PERIODIC_X) || defined(PERIODIC_Y)
		make_periodic(&new_position);
		#endif
		double new_phi=angle_twopi(particles[i].phi+max_rotation*(dsfmt_genrand_open_close(&rng)-0.5));

		int du_ext=0;
		int du_int=0;
		if(mc_energy_change(new_position,new_phi,i,&du_ext,&du_int) && ( du_ext*ENERGY_WELL_DEPTH+du_int*ENERGY_BOND <= 0 || dsfmt_genrand_open_close(&rng) <= mc_acceptance_probability(du_ext, du_int) ) ){
			particles[i].position=new_position;
			particles[i].phi=new_phi;
		}
	}
	return 1.0;
}

double mc_run(int steps){
	double md=0.1;
	double mr=M_PI/10.0;
	log_enqueue(0,false);
	for(int i=0;i<steps;){
		for(int j=0;j<LOGGING_INTERVAL && i<steps;++j){
			monte_carlo_step(md,mr);
			++i;
		}
		log_enqueue(i,i==steps);
	}
	return 1.0;
}

bool mc_energy_change(vector2d new_position, double new_phi, int i, int *e_ext, int *e_int){
	*e_ext=0;
	*e_int=0;
	double d=0;
	for(int j=0;j<NUMBER_OF_PARTICLES;++j){
		if(j==i) continue;
		distance(particles[i].position,particles[j].position,&d);
		if(d<COLLOID_DIAMETER){
			return false;
		}
	}
	return true;
}

double mc_acceptance_probability(int du_ext, int du_int){
	return acceptance_probabilities_well[1+du_ext]*acceptance_probabilities_bonds[3+du_int];
}

void mc_init_acceptance_probabilities(double kbt){
	for(int du_ext=-1;du_ext<=1;++du_ext){
		acceptance_probabilities_well[1+du_ext]=exp(-du_ext*ENERGY_WELL_DEPTH/kbt);
	}
	for(int du_int=-3;du_int<=3;++du_int){
		acceptance_probabilities_bonds[3+du_int]=exp(-du_int*ENERGY_BOND/kbt);
	}
}

void mc_init_particles(void){
	double d;
	bool collision=false;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		do{
			particles[i].position[0]=SIZE_X*dsfmt_genrand_open_close(&rng);
			particles[i].position[1]=SIZE_Y*dsfmt_genrand_open_close(&rng);
			particles[i].phi=2.0*M_PI*dsfmt_genrand_open_close(&rng);
			collision=false;
			for(int j=0;j<i;++j){
				distance(particles[i].position,particles[j].position,&d);
				if(d<COLLOID_DIAMETER){
					collision=true;
					break;
				}
			}
		}while(collision==true);
	}
}

void mc_init(double kbt){
	mc_init_particles();
	mc_init_acceptance_probabilities(kbt);
}

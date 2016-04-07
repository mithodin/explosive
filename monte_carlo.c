#include <math.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "dSFMT/dSFMT.h"
#include "monte_carlo.h"
#include "colloid.h"
#include "parameters.h"
#include "distance.h"

double monte_carlo_step(double, double);
bool mc_energy_change(vector2d, double, int, int*, int*);
double mc_acceptance_probability(int, int);
void mc_init_particles(void);
void mc_init_acceptance_probabilities(double);

double acceptance_probabilities_bonds[7];
double acceptance_probabilities_well[3];
Colloid *particles;

double monte_carlo_step(double max_displacement, double max_rotation){
	for(int i=1;i<num_particles;++i){
		vector2d new_position;
		new_position[0]=particles[i].position[0]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		new_position[1]=particles[i].position[1]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		double new_phi=angle_twopi(particles[i].phi+max_rotation*(dsfmt_genrand_open_close(&rng)-0.5));

		int du_ext=0;
		int du_int=0;
		if(mc_energy_change(new_position,new_phi,i,&du_ext,&du_int) && ( du_ext*energy_well_depth+du_int*energy_bond <= 0 || dsfmt_genrand_open_close(&rng) <= mc_acceptance_probability(du_ext, du_int) ) ){
			particles[i].position=new_position;
			particles[i].phi=new_phi;
		}
	}
	return 1.0;
}

double mc_run(int steps){
	double md=0.01;
	double mr=M_PI/10.0;
	for(int i=0;i<steps;++i){
		monte_carlo_step(md,mr);
	}
	return 1.0;
}

bool mc_energy_change(vector2d new_position, double new_phi, int i, int *e_ext, int *e_int){
	if(new_position[0] > SIZE_X || new_position[0] < 0 || new_position[1] > SIZE_Y || new_position[1] < 0){
		return false;
	}
	*e_ext=0;
	*e_int=0;
	double d=0;
	distance(particles[0].position,particles[1].position,&d);
	return d>colloid_diameter;
}

double mc_acceptance_probability(int du_ext, int du_int){
	return acceptance_probabilities_well[1+du_ext]*acceptance_probabilities_bonds[3+du_int];
}

void mc_init_acceptance_probabilities(double kbt){
	for(int du_ext=-1;du_ext<=1;++du_ext){
		acceptance_probabilities_well[1+du_ext]=exp(-du_ext*energy_well_depth/kbt);
	}
	for(int du_int=-3;du_int<=3;++du_int){
		acceptance_probabilities_bonds[3+du_int]=exp(-du_int*energy_bond/kbt);
	}
}

void mc_init_particles(void){
	particles=calloc(2,sizeof(Colloid));
	particles[0].position[0]=0.0;
	particles[0].position[1]=0.0;
	particles[0].phi=0.0;
	particles[1].position[0]=1.0;
	particles[1].position[1]=1.0;
	particles[1].phi=M_PI;
}

void mc_init(double kbt){
	mc_init_particles();
	mc_init_acceptance_probabilities(kbt);
}

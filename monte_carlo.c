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
#include "substrate.h"

double monte_carlo_step(double, double);
bool mc_energy_change(Colloid *, int);
double mc_acceptance_probability(int, int);
void mc_init_particles(void);
void mc_init_acceptance_probabilities(double);

double acceptance_probabilities_bonds[7];
double acceptance_probabilities_well[3];

double monte_carlo_step(double max_displacement, double max_rotation){
	int accept=0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		Colloid new=EMPTY_COLLOID;
		new.position[0]=particles[i].position[0]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		#ifndef PERIODIC_X
		if(new.position[0] < 0 || new.position[0] > SIZE_X){
			continue;
		}
		#endif
		new.position[1]=particles[i].position[1]+max_displacement*(dsfmt_genrand_open_close(&rng)-0.5);
		#ifndef PERIODIC_Y
		if(new.position[1] < 0 || new.position[1] > SIZE_Y){
			continue;
		}
		#endif
		#if defined(PERIODIC_X) || defined(PERIODIC_Y)
		make_periodic(&(new.position));
		#endif
		new.phi=angle_twopi(particles[i].phi+max_rotation*(dsfmt_genrand_open_close(&rng)-0.5));

		if(mc_energy_change(&new,i) && ( (new.external_energy-particles[i].external_energy)*ENERGY_WELL_DEPTH+(new.internal_energy-particles[i].internal_energy)*ENERGY_BOND <= 0 || dsfmt_genrand_open_close(&rng) <= mc_acceptance_probability((new.external_energy-particles[i].external_energy), (new.internal_energy-particles[i].internal_energy)) ) ){
			particles[i].position=new.position;
			particles[i].phi=new.phi;
			particles[i].external_energy=new.external_energy;
			particles[i].internal_energy=new.internal_energy;
			particles[i].below->above=particles[i].above;
			particles[i].above->below=particles[i].below;
			insert_sorted_y(&(particles[i]),particles[i].below);

			//remove old bonding partners
			for(int k=0;k<3;++k){
				Colloid *partner=particles[i].bonding_partner[k];
				if(partner!=NULL){
					partner->bonding_partner[particles[i].bond_site[k]]=NULL;
					partner->internal_energy++;
				}
			}
			//make new bonds
			for(int k=0;k<3;++k){
				Colloid *partner=new.bonding_partner[k];
				particles[i].bonding_partner[k]=partner;
				particles[i].bond_site[k]=new.bond_site[k];
				if(partner!=NULL){
					partner->bonding_partner[new.bond_site[k]]=&(particles[i]);
					partner->bond_site[new.bond_site[k]]=k;
					partner->internal_energy--;
				}
			}
			++accept;
		}
	}
	return 1.0*accept/NUMBER_OF_PARTICLES;
}

double mc_run(int steps){
	char percent_complete[103];
	double md=0.1;
	double mr=M_PI/10.0;
	double complete=0.0;
	log_enqueue(0,false);
	for(int i=0;i<steps;){
		for(int j=0;j<LOGGING_INTERVAL && i<steps;++j){
			monte_carlo_step(md,mr);
			++i;
		}
		complete=1.0*i/steps;
		mkpercent(percent_complete,103,complete);
		printf("\r> running %s %3d%% complete",percent_complete,(int)floor(100*complete));
		fflush(NULL);
		log_enqueue(i,i==steps);
	}
	printf("\n");
	return 1.0;
}

bool mc_energy_change(Colloid *new, int i){
	new->external_energy=external_energy(new->position);
	new->internal_energy=0;
	bool collision;
	int site1;
	int site2;
	Colloid *current=particles[i].above;
	while(distance_y(new->position[1],current->position[1]) <= COLLOID_MIN_BONDING_DISTANCE){
		if(colloid_bonded(current,new,&collision,&site1,&site2) && !collision){
			new->internal_energy--;
			new->bonding_partner[site2]=current;
			new->bond_site[site2]=site1;
		}else if(collision){
			return false;
		}
		current=current->above;
	}
	current=particles[i].below;
	while(distance_y(current->position[1],new->position[1]) <= COLLOID_MIN_BONDING_DISTANCE){
		if(colloid_bonded(current,new,&collision,&site1,&site2) && !collision){
			new->internal_energy--;
			new->bonding_partner[site2]=current;
			new->bond_site[site2]=site1;
		}else if(collision){
			return false;
		}
		current=current->below;
	}
	//naive implementation
	/*for(int j=0;j<NUMBER_OF_PARTICLES;++j){
		if(j==i) continue;
		distance(new_position,particles[j].position,&d);
		if(d<COLLOID_DIAMETER){
			return false;
		}
	}*/
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
	init_ysorted_list();
	init_bonding_partners();
}

void mc_init(double kbt){
	mc_init_particles();
	mc_init_acceptance_probabilities(kbt);
}

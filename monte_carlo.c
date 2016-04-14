/** @file */

#include <math.h>
#include <sys/time.h>
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
#include "clusters.h"

double monte_carlo_step(void);
bool mc_energy_change(Colloid *, int);
double mc_acceptance_probability(int, int);
void mc_init_particles(void);
void mc_init_acceptance_probabilities(double);
void mc_init_max_displacement(double);
double timediff_seconds(struct timeval *, struct timeval *);

double acceptance_probabilities_bonds[7]; /**< pre-calculated acceptance probabilities for breaking and making bonds */
double acceptance_probabilities_well[3]; /**< pre-calculated acceptance probabilities for entering and leaving a well on the substrate */
double max_displacement=0.1; /**< maximum displacement in either coordinate during one mc move */
double max_rotation=M_PI; /**< maximum rotation during one mc move */ 

/**
 * Do one monte carlo step (try one move for every particle)
 *
 * @return The acceptance probability
 */
double monte_carlo_step(void){
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

/**
 * Run some monte carlo steps
 *
 * @param steps How many steps to run
 * @param log Whether or not to log this run in the results file
 * @return The acceptance probability
 */
double mc_run(int steps, bool log){
	char percent_complete[103];
	double complete=0.0;
	double acceptance_probability=0.0;
	unsigned long runtime;
	int largest_cluster;
	struct timeval t0,t1;
	if(log){
		gettimeofday(&t0,NULL);
		log_enqueue(0,false,(unsigned long)(t0.tv_sec-sim_start_time.tv_sec),0);
	}
	for(int i=0;i<steps;){
		for(int j=0;j<LOGGING_INTERVAL && i<steps;++j){
			acceptance_probability+=monte_carlo_step();
			++i;
		}
		if(log){
			gettimeofday(&t1,NULL);
			runtime=(unsigned long)(t1.tv_sec-sim_start_time.tv_sec);
			largest_cluster=largest_cluster_size();
			log_enqueue(i,i==steps,runtime,largest_cluster);
			complete=1.0*i/steps;
			mkpercent(percent_complete,103,complete);
			printf("\r> running %lds %s %3d%% complete. ETA: %.0fs       ",runtime,percent_complete,(int)floor(100*complete),1.0*(steps-i)/LOGGING_INTERVAL*timediff_seconds(&t1,&t0));
			fflush(NULL);
			t0=t1;
		}
	}
	if(log) printf("\n");
	return acceptance_probability/steps;
}

/**
 * Calculate the internal and external energy for the new particle
 *
 * @param new The new colloid
 * @param i Index of the soon-to-be-updated particle
 * @return Was there a collision?
 */
bool mc_energy_change(Colloid *new, int i){
	new->external_energy=external_energy(new->position);
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

/**
 * Calculate the acceptance probability for a given energy change
 *
 * @param du_ext Has the particle entered (-1) or left (+1) a well on the substrate?
 * @param du_int Have the bonds of the particle changed?
 * @return a number between 0.0 and 1.0
 */
double mc_acceptance_probability(int du_ext, int du_int){
	return acceptance_probabilities_well[1+du_ext]*acceptance_probabilities_bonds[3+du_int];
}

/**
 * Precalculate the acceptance probabilities so we don't have to evaluate millions of exp()s.
 *
 * @param kbt The temperature of the system
 */
void mc_init_acceptance_probabilities(double kbt){
	for(int du_ext=-1;du_ext<=1;++du_ext){
		acceptance_probabilities_well[1+du_ext]=exp(-du_ext*ENERGY_WELL_DEPTH/kbt);
	}
	for(int du_int=-3;du_int<=3;++du_int){
		acceptance_probabilities_bonds[3+du_int]=exp(-du_int*ENERGY_BOND/kbt);
	}
}

/**
 * Place particles randomly in the box
 */
void mc_init_particles(void){
	printf("> initializing particles... ");
	fflush(NULL);
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
		}while(collision);
		particles[i].external_energy=external_energy(particles[i].position);
		particles[i].particles_index=i;
	}
	init_ysorted_list();
	init_bonding_partners();
	printf("done.\n");
}

/**
 * Tune maximum displacement and rotation so that a target acceptance rate is met.
 *
 * @param target_acceptance_rate The required acceptance rate
 */
void mc_init_max_displacement(double target_acceptance_rate){
	printf("> initializing maximum displacement... ");
	fflush(NULL);
	double md_tmp=max_displacement;
	double tar_sqrt=sqrt(target_acceptance_rate);
	double acceptance_rate=mc_run(100,false);
	int i=0,j=0;
	while( j < 10 && fabs(acceptance_rate-target_acceptance_rate)/target_acceptance_rate > 0.01){
		md_tmp=max_displacement;
		max_displacement=0.0;
		acceptance_rate=mc_run(100,false);
		i=0;
		while(max_rotation<2.0*M_PI && i < 100 && fabs(acceptance_rate-tar_sqrt)/tar_sqrt > 0.01){
			max_rotation*=acceptance_rate/tar_sqrt;
			acceptance_rate=mc_run(100,false);
			++i;
		}
		max_displacement=md_tmp;
		acceptance_rate=mc_run(100,false);
		i=0;
		while( i < 100 && fabs(acceptance_rate-target_acceptance_rate)/target_acceptance_rate > 0.01){
			max_displacement*=acceptance_rate/target_acceptance_rate;
			acceptance_rate=mc_run(100,false);
			++i;
		}
		++j;
	}
	if(j==10 || i==100){
		printf("[too many iterations] ");
	}
	printf("done.\n");
}

/**
 * Initialize the monte carlo subsystem including:
 *  - the substrate
 *  - the particles
 *  - the acceptance probabilities
 *  - the maximum displacement and rotation
 *
 * @param kbt The temperature of the system
 */
void mc_init(double kbt){
	init_substrate();
	mc_init_particles();

	//thermalize
	mc_init_acceptance_probabilities(2.0);
	mc_run(100,false);

	mc_init_acceptance_probabilities(kbt);
	mc_init_max_displacement(0.5);
}

/**
 * Calculate the difference between to timevals
 *
 * @param t1 The end time
 * @param t0 The start time
 * @return Difference between t0 and t1 in seconds
 */
double timediff_seconds(struct timeval *t1, struct timeval *t0){
	return 1.0*(t1->tv_sec-t0->tv_sec)+1e-6*(t1->tv_usec-t0->tv_usec);
}

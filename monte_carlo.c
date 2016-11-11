/** @file */
#include <math.h>
#include <sys/time.h>
#include <x86intrin.h>
#include <stdbool.h>
#include <string.h>
#include <stdlib.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "geometry.h"
#include "colloid.h"
#include "logger.h"
#include "globals.h"
#include "monte_carlo.h"
#include "substrate.h"
#include "clusters.h"

#ifdef CONTINUE
#include <hdf5.h>
#include <hdf5_hl.h>

hid_t configuration;
hid_t conf_group;
unsigned int old_number_of_particles;
void *old_positions;
bool oldpos_float;
  #ifndef OLD_LOGFILE
extern hid_t logfile;
  #endif
#endif

double monte_carlo_step(void);
bool mc_energy_change(Colloid *, int);
double mc_acceptance_probability(double, int);
bool mc_init_particles(void);
void mc_init_acceptance_probabilities(double);
bool mc_init_max_displacement(double);
double timediff_seconds(struct timeval *, struct timeval *);
unsigned long time_hm(unsigned long, unsigned long *, unsigned long *);
void shuffle_float(size_t, size_t, float *);
void shuffle_double(size_t, size_t, double *);

double acceptance_probabilities_bonds[7]; /**< pre-calculated acceptance probabilities for breaking and making bonds */
double temperature; /**< System temperature in kbT */
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
		new.position.v=_mm_add_pd(particles[i].position.v,_mm_set_pd(max_displacement*(dsfmt_genrand_open_close(&rng)-0.5),max_displacement*(dsfmt_genrand_open_close(&rng)-0.5)));
		#ifndef PERIODIC_X
		if(new.position.c.x < 0 || new.position.c.x > SIZE_X){
			continue;
		}
		#endif
		#ifndef PERIODIC_Y
		if(new.position.c.y < 0 || new.position.c.y > SIZE_Y){
			continue;
		}
		#endif
		#if defined(PERIODIC_X) || defined(PERIODIC_Y)
		new.position=make_periodic(new.position);
		#endif
		new.phi=angle_twopi(particles[i].phi+max_rotation*(dsfmt_genrand_open_close(&rng)-0.5));

		if(mc_energy_change(&new,i) && ( (new.external_energy-particles[i].external_energy)+(new.internal_energy-particles[i].internal_energy)*ENERGY_BOND <= 0 || dsfmt_genrand_open_close(&rng) <= mc_acceptance_probability((new.external_energy-particles[i].external_energy), (new.internal_energy-particles[i].internal_energy)) ) ){
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
	double acceptance_probability=0.0,bonding_probability;
	unsigned long runtime,rth,rtm,hours,minutes,seconds;
	int largest_cluster;
	struct timeval t0,t1;
	if(log){
		gettimeofday(&t0,NULL);
		log_enqueue(0,false,(unsigned long)(t0.tv_sec-sim_start_time.tv_sec),0,0);
	}
	for(int i=0;i<steps;){
		for(int j=0;j<LOGGING_INTERVAL && i<steps;++j){
			acceptance_probability+=monte_carlo_step();
			++i;
		}
		if(log){
			gettimeofday(&t1,NULL);
			runtime=time_hm((unsigned long)(t1.tv_sec-sim_start_time.tv_sec),&rtm,&rth);
			largest_cluster=largest_cluster_size(&bonding_probability);
			log_enqueue(i,i==steps,runtime,largest_cluster,bonding_probability);
			complete=1.0*i/steps;
			mkpercent(percent_complete,103,complete);
			seconds=time_hm((unsigned long)(1.0*(steps-i)/LOGGING_INTERVAL*timediff_seconds(&t1,&t0)),&minutes,&hours);
			printf("\r> running %ldh %2ldm %2lds %s %3d%% complete. ETA: %2ldh %2ldm %2lds       ",rth,rtm,runtime,percent_complete,(int)floor(100*complete),hours,minutes,seconds);
			fflush(NULL);
			t0=t1;
		}
	}
	if(log) printf("\n");
	return acceptance_probability/steps;
}

/**
 * Run monte carlo steps until the system is equilibrated (see config.h)
 */
void mc_equilibrate(void){
	FILE *energy_log_debug=fopen("energy_dbg.dat","w");
	double energy_buffer[EQUILIBRATION_SMOOTHING_STEP]={0};
	double energy_buffer_smoothed[EQUILIBRATION_SMOOTHING_STEP]={0};
	double energy_cumulative=0;
	int energy_buffer_in=0;
	int count_loops=0;
	bool first_loop_done=false;
	bool equilibrated=false;
	printf("> Begin equilibration...\n");
	while(!equilibrated){
		mc_run(EQUILIBRATION_INNER_LOOP,false);
		count_loops++;
		energy_cumulative-=energy_buffer[energy_buffer_in]/EQUILIBRATION_SMOOTHING_STEP;
		energy_buffer[energy_buffer_in]=0;
		for(int i=0;i<NUMBER_OF_PARTICLES;++i){
			energy_buffer[energy_buffer_in]+=particles[i].external_energy+particles[i].internal_energy*ENERGY_BOND/2.0;
		}
		energy_cumulative+=energy_buffer[energy_buffer_in]/EQUILIBRATION_SMOOTHING_STEP;
		double m=2.0*EQUILIBRATION_THRESHOLD_SLOPE; //definitely above the threshold
		if(!first_loop_done){
			energy_buffer_smoothed[energy_buffer_in]=0;
			for(int i=0;i<=energy_buffer_in;++i){
				energy_buffer_smoothed[energy_buffer_in]+=energy_buffer[i];
			}
			energy_buffer_smoothed[energy_buffer_in]/=energy_buffer_in+1;
			if(energy_buffer_in+1==EQUILIBRATION_SMOOTHING_STEP){
				first_loop_done=true;
			}
			if(energy_buffer_in > 0){
				m=(energy_buffer_smoothed[0]-energy_buffer_smoothed[energy_buffer_in])/(EQUILIBRATION_INNER_LOOP*energy_buffer_in);
			}
		}
		else if(first_loop_done){
			m=fabs(energy_buffer_smoothed[energy_buffer_in]-energy_cumulative)/(EQUILIBRATION_SMOOTHING_STEP*EQUILIBRATION_INNER_LOOP);
			energy_buffer_smoothed[energy_buffer_in]=energy_cumulative;
		}
		if( m < EQUILIBRATION_THRESHOLD_SLOPE ){
			printf("\r> system equilibrated after %1.5e MC steps                                                    \n",1.0*count_loops*EQUILIBRATION_INNER_LOOP);
			equilibrated=true;
		}
		else{
			printf("\r> equilibrating (current slope %1.7e, threshold %1.2e) [%1.5e MC steps]          ",m,EQUILIBRATION_THRESHOLD_SLOPE,1.0*count_loops*EQUILIBRATION_INNER_LOOP);
			fflush(NULL);
		}
		fprintf(energy_log_debug,"%d\t%1.5e\t%1.5e\t%1.5e\n",count_loops*EQUILIBRATION_INNER_LOOP,energy_buffer[energy_buffer_in],energy_buffer_smoothed[energy_buffer_in],m);
		fflush(energy_log_debug);
		energy_buffer_in=(energy_buffer_in+1)%EQUILIBRATION_SMOOTHING_STEP;
	}
	fclose(energy_log_debug);

	double runtime=0.0,bonding_probability;
	double largest_cluster=largest_cluster_size(&bonding_probability);
	log_enqueue(0,false,runtime,largest_cluster,bonding_probability);

	log_checkpoint(); //if it fails for some reason, it's a problem, but not a fatal one.
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
	while(distance_y(new->position.c.y,current->position.c.y) <= COLLOID_MIN_BONDING_DISTANCE){
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
	while(distance_y(current->position.c.y,new->position.c.y) <= COLLOID_MIN_BONDING_DISTANCE){
		if(colloid_bonded(current,new,&collision,&site1,&site2) && !collision){
			new->internal_energy--;
			new->bonding_partner[site2]=current;
			new->bond_site[site2]=site1;
		}else if(collision){
			return false;
		}
		current=current->below;
	}
	return true;
}

/**
 * Calculate the acceptance probability for a given energy change
 *
 * @param du_ext Has the particle entered (-1) or left (+1) a well on the substrate?
 * @param du_int Have the bonds of the particle changed?
 * @return a number between 0.0 and 1.0
 */
double mc_acceptance_probability(double du_ext, int du_int){
	return exp(-du_ext/temperature)*acceptance_probabilities_bonds[3+du_int];
}

/**
 * Precalculate the acceptance probabilities so we don't have to evaluate millions of exp()s.
 *
 * @param kbt The temperature of the system
 */
void mc_init_acceptance_probabilities(double kbt){
	temperature=kbt;
	for(int du_int=-3;du_int<=3;++du_int){
		acceptance_probabilities_bonds[3+du_int]=exp(-du_int*ENERGY_BOND/kbt);
	}
}

/**
 * Place particles randomly in the box
 *
 * @return Were the particles successfully initialized?
 */
bool mc_init_particles(void){
	if(NUMBER_OF_PARTICLES*M_PI*COLLOID_DIAMETER/4.0 > SIZE_X*SIZE_Y){
		printf("> error: too many particles\n");
		return false;
	}
	printf("> initializing particles... ");
	fflush(NULL);
#ifdef CONTINUE
	if(oldpos_float){
		shuffle_float((size_t)old_number_of_particles,3,(float*)old_positions);
		for(int i=0;i<old_number_of_particles && i < NUMBER_OF_PARTICLES;++i){
			particles[i].position.c.x=((float *)old_positions)[3*i];
			particles[i].position.c.y=((float *)old_positions)[3*i+1];
			particles[i].phi=((float *)old_positions)[3*i+1];
			particles[i].external_energy=external_energy(particles[i].position);
			particles[i].particles_index=i;
		}
	}else{
		shuffle_double((size_t)old_number_of_particles,3,(double*)old_positions);
		for(int i=0;i<old_number_of_particles && i < NUMBER_OF_PARTICLES;++i){
			particles[i].position.c.x=((double *)old_positions)[3*i];
			particles[i].position.c.y=((double *)old_positions)[3*i+1];
			particles[i].phi=((double *)old_positions)[3*i+1];
			particles[i].external_energy=external_energy(particles[i].position);
			particles[i].particles_index=i;
		}
	}
	free(old_positions);
#elif PARTICLES_INIT_RANDOM == 0
	int particles_num_cols,particles_num_rows;
	particles_num_cols=(int)ceil(sqrt(NUMBER_OF_PARTICLES*SIZE_X*1.0/SIZE_Y));
	particles_num_rows=(int)ceil(1.0*NUMBER_OF_PARTICLES/particles_num_cols);
	double step_x=SIZE_X/particles_num_cols;
	double step_y=SIZE_Y/particles_num_rows;
#endif
	double d;
	bool collision=false;
#ifdef CONTINUE
	for(int i=old_number_of_particles;i<NUMBER_OF_PARTICLES;++i){
#else
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
#endif
#if PARTICLES_INIT_RANDOM == 1 || defined(CONTINUE)
		do{
			particles[i].position.v=_mm_set_pd(SIZE_Y*dsfmt_genrand_open_close(&rng),SIZE_X*dsfmt_genrand_open_close(&rng));
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
#else //init in a square lattice
		particles[i].position.v=_mm_set_pd(((i/particles_num_cols)+0.5)*step_y,((i%particles_num_cols)+0.5)*step_x);
		particles[i].phi=2.0*M_PI*dsfmt_genrand_open_close(&rng);
		for(int j=0;j<i;++j){
			distance(particles[i].position,particles[j].position,&d);
			if(d<COLLOID_DIAMETER){
				printf("> Error: Too many particles\n");
				return false;
			}
		}
#endif
		particles[i].external_energy=external_energy(particles[i].position);
		particles[i].particles_index=i;
	}
	init_ysorted_list();
	init_bonding_partners();
	printf("done.\n");
	return true;
}

/**
 * Tune maximum displacement and rotation so that a target acceptance rate is met.
 *
 * @param target_acceptance_rate The required acceptance rate
 */
bool mc_init_max_displacement(double target_acceptance_rate){
	printf("> initializing maximum displacement... ");
	fflush(NULL);
	double md_tmp=max_displacement;
	double tar_sqrt=sqrt(target_acceptance_rate);
	double acceptance_rate=mc_run(100,false);
	int i=0,j=0;
	printf("dmax: %8.5f amax: %7.5f  ",max_displacement,max_rotation);
	fflush(NULL);
	/*
	while( j < 10 && fabs(acceptance_rate-target_acceptance_rate)/target_acceptance_rate > 0.01){
		md_tmp=max_displacement;
		max_displacement=0.0;
		acceptance_rate=mc_run(100,false);
		i=0;
		while(i < 100 && fabs(acceptance_rate-tar_sqrt)/tar_sqrt > 0.01){
			max_rotation*=acceptance_rate/tar_sqrt;
			if( max_rotation > 2.0*M_PI ){
				max_rotation = 2.0*M_PI;
				break;
			}
			printf("\r> initializing maximum displacement... dmax: %8.5f amax: %7.5f ar: %3.0f%% ",md_tmp,max_rotation,acceptance_rate*acceptance_rate*100);
			fflush(NULL);
			acceptance_rate=mc_run(100,false);
			++i;
		}
		max_displacement=md_tmp;
		acceptance_rate=mc_run(100,false);
		printf("\r> initializing maximum displacement... dmax: %8.5f amax: %7.5f ar: %3.0f%% ",max_displacement,max_rotation,acceptance_rate*100);
		fflush(NULL);
		i=0;
		while( i < 100 && fabs(acceptance_rate-target_acceptance_rate)/target_acceptance_rate > 0.01){
			max_displacement*=acceptance_rate/target_acceptance_rate;
			if( max_displacement > SIZE_X ){
				max_displacement = SIZE_X;
				break;
			}
			printf("\r> initializing maximum displacement... dmax: %8.5f amax: %7.5f ar: %3.0f%% ",max_displacement,max_rotation,acceptance_rate*100);
			fflush(NULL);
			acceptance_rate=mc_run(100,false);
			++i;
		}
		++j;
	}
	if(j==10 || i==100){
		printf("[too many iterations] ");
	}
	*/
	while( fabs(acceptance_rate-target_acceptance_rate)/target_acceptance_rate > 0.01 ){
		acceptance_rate = mc_run(100,false);
		md_tmp = max_displacement*acceptance_rate/target_acceptance_rate;
		max_displacement = 0.0;
		acceptance_rate = mc_run(200,false);
		max_rotation *= acceptance_rate/tar_sqrt;
		if( max_rotation > 2.0*M_PI ){
			max_rotation = 2.0*M_PI;
		}
		max_displacement = md_tmp;
		acceptance_rate = mc_run(100,false);
		max_displacement *= acceptance_rate/target_acceptance_rate;
		if( max_displacement > SIZE_X ){
			max_displacement = SIZE_X;
		}
		printf("\r> initializing maximum displacement... dmax: %8.5f amax: %7.5f ar: %3.0f%% ",max_displacement,max_rotation,acceptance_rate*100);
		fflush(NULL);
		if( j == 1000 ){ //this is going nowhere. starting over.
			max_displacement=0.1;
			max_rotation=M_PI; 
			j = 0;
			printf("\n> starting over.\n");
		}
	}
	printf("done.\n");
	return log_max_displacement();
}

/**
 * Initialize the monte carlo subsystem including:
 *  - the substrate
 *  - the particles
 *  - the acceptance probabilities
 *  - the maximum displacement and rotation
 *
 * @param kbt The temperature of the system
 * @return Was the monte carlo subsystem successfully intialized?
 */
bool mc_init(double kbt){
#ifdef CONTINUE
	herr_t status;
	hsize_t number_of_fields,number_of_records;
  #ifdef OLD_LOGFILE
	configuration = H5Fopen(OLD_LOGFILE, H5F_ACC_RDONLY, H5P_DEFAULT);
  #else
	configuration = logfile;
  #endif
  	if( configuration < 0 ){
		printf("> Could not open old configuration file\n");
		return false;
	}

	conf_group = H5Gopen(configuration, OLD_LOGFILE_GROUP, H5P_DEFAULT);
	if( conf_group < 0 ){
		printf("> Could not open group.\n");
		return false;
	}
#endif
	if( !init_substrate() ){
		printf("> init_substrate failed\n");
		return false;
	};
#ifdef CONTINUE
	status = H5LTget_attribute_uint(conf_group,OLD_LOGFILE_GROUP,"number-of-particles",&old_number_of_particles);
	if( status < 0 ){ printf("> Error reading number of particles from old file\n"); return false; }
	status = H5LTget_attribute_double(conf_group,OLD_LOGFILE_GROUP,"max-displacement",&max_displacement);
	if( status < 0 ){ printf("> Error reading max displacement from old file\n"); return false; }
	status = H5LTget_attribute_double(conf_group,OLD_LOGFILE_GROUP,"max-rotation",&max_rotation);
	if( status < 0 ){ printf("> Error reading max rotation from old file\n"); return false; }

	status = H5TBget_table_info(conf_group,"simulation_frames",&number_of_fields,&number_of_records);
	if( status < 0 ){ printf("> Error reading table info from old file\n"); return false; }

	char **field_names = calloc(number_of_fields,sizeof(char *));
	for(int i=0;i<number_of_fields;++i){
		field_names[i]=calloc(50,sizeof(char));
	}
	size_t *field_sizes = calloc(number_of_fields,sizeof(size_t));
	size_t *field_offsets = calloc(number_of_fields,sizeof(size_t));
	size_t type_size=0;

	status = H5TBget_field_info(conf_group,"simulation_frames",field_names,field_sizes,field_offsets,&type_size);
	if( status < 0){
		printf("> Error reading field info.\n");
		return false;
	}

	size_t position=0;
	for(int i=0;i<number_of_fields;++i){
		if(strcmp(field_names[i],"position")==0){ position=field_offsets[i]; }
	}

	void *record_buffer=malloc(type_size);
	status = H5TBread_records(conf_group, "simulation_frames", number_of_records-1, 1, type_size, field_offsets, field_sizes, record_buffer);
	if( status < 0 ){ printf("> Error reading old configuration\n"); return false; }

	for(int i=0;i<number_of_fields;++i){
		if(strcmp(field_names[i],"position")==0){
			oldpos_float=field_sizes[i]/old_number_of_particles/3==sizeof(float);
		}
	}

	old_positions = calloc(old_number_of_particles, oldpos_float?3*sizeof(float):3*sizeof(double));
	memcpy(old_positions,record_buffer+position,old_number_of_particles*3*(oldpos_float?sizeof(float):sizeof(double)));

	free(record_buffer);
	status = H5Gclose(conf_group);
	if(status < 0){ printf("> Error closing the oldconf group\n"); return false; }

  #ifdef OLD_LOGFILE
	status = H5Fclose(configuration);
	if(status < 0){ printf("> Error closing the oldconf file\n"); return false; }
  #endif
#endif
	if( !mc_init_particles() ){
		return false;
	}

#ifndef CONTINUE
	//thermalize
	mc_init_acceptance_probabilities(100.0);
	mc_run(1000,false);
#endif

	mc_init_acceptance_probabilities(kbt);
	return mc_init_max_displacement(0.5);
}

//Helper functions

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

/**
 * Calculate hours and minutes from seconds
 */
unsigned long time_hm(unsigned long seconds, unsigned long *minutes, unsigned long *hours){
	*hours=seconds/3600;
	seconds-=*hours*3600;
	*minutes=seconds/60;
	seconds-=*minutes*60;
	return seconds;
}

/**
 * Shuffle a size*elem-length array of floats. Blocks of length elem will be kept intact
 *
 * @param size The number of elem-length blocks of floats in the array
 * @param elem The length of the elementary blocks (will be swapped as one)
 */
void shuffle_float(size_t size, size_t elem, float *array){
	float *tmp=calloc(elem,sizeof(float));
	for(int i=(int)size-1;i>=0;--i){
		int random_index = (int)floor(dsfmt_genrand_close_open(&rng)*i);
		if(random_index != i){
			memcpy(tmp,array+random_index*elem,elem*sizeof(float));
			memcpy(array+random_index*elem,array+i*elem,elem*sizeof(float));
			memcpy(array+i*elem,tmp,elem*sizeof(float));
		}
	}
	free(tmp);
}

/**
 * Shuffle a size*elem-length array of doubles. Blocks of length elem will be kept intact
 *
 * @param size The number of elem-length blocks of doubles in the array
 * @param elem The length of the elementary blocks (will be swapped as one)
 */
void shuffle_double(size_t size, size_t elem, double *array){
	double *tmp=calloc(elem,sizeof(double));
	for(int i=(int)size-1;i>=0;--i){
		int random_index = (int)floor(dsfmt_genrand_close_open(&rng)*i);
		if(random_index != i){
			memcpy(tmp,array+random_index*elem,elem*sizeof(double));
			memcpy(array+random_index*elem,array+i*elem,elem*sizeof(double));
			memcpy(array+i*elem,tmp,elem*sizeof(double));
		}
	}
	free(tmp);
}

/** @file */

#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <x86intrin.h>
#include <stdbool.h>
#include <unistd.h>
#include <stdlib.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "geometry.h"
#include "colloid.h"
#include "simulation_frame.h"
#include "globals.h"

hid_t logfile;
hid_t group;
char directory_name[102];
char simulation_frames_location[120];
size_t simframe_size;
size_t simframe_offsets[7];

char cluster_size_distribution_location[120];
size_t cluster_bin_size;
size_t cluster_bin_offsets[3];

char patch_storage_location[120];

/**
 * Initialize the hdf5 log file. Set up groups and tables as well as attributes
 * @return Did the initialization succeed?
 */
bool h5log_init(void){
	herr_t status;

	if( access( LOGFILE, F_OK ) == 0 ){
		logfile = H5Fopen(LOGFILE, H5F_ACC_RDWR, H5P_DEFAULT);
		if(logfile < 0){
			return false;
		}
	}else{ //file does not exist
		logfile = H5Fcreate(LOGFILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if(logfile < 0){
			printf("> H5Log cannot open a file to log this simulation\n");
			return false;
		}
	}

	//create group
	snprintf(directory_name, sizeof(directory_name), "/%.100s", SIMULATION_SHORT_NAME);

	int dir_exists_index=0;
	while(H5Lexists(logfile, directory_name, H5P_DEFAULT)){
		snprintf(directory_name, sizeof(directory_name), "/%.94s_%05d", SIMULATION_SHORT_NAME, dir_exists_index++);
	}
	snprintf(simulation_frames_location, sizeof(simulation_frames_location), "%.101s/simulation_frames", directory_name);
	snprintf(cluster_size_distribution_location, sizeof(cluster_size_distribution_location), "%.101s/cluster_size", directory_name);
	snprintf(patch_storage_location, sizeof(patch_storage_location), "%.101s/random_patches", directory_name);
	group = H5Gcreate(logfile, directory_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//create group attributes (stores simulation parameters)
	status = H5LTset_attribute_string(group, directory_name, "simulation-name", SIMULATION_NAME);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const unsigned int number_of_particles=NUMBER_OF_PARTICLES;
	status = H5LTset_attribute_uint(group, directory_name, "number-of-particles", &number_of_particles, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double size_x=SIZE_X;
	status = H5LTset_attribute_double(group, directory_name, "size-x", &size_x, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double size_y=SIZE_Y;
	status = H5LTset_attribute_double(group, directory_name, "size-y", &size_y, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	#ifdef PERIODIC_X
	status = H5LTset_attribute_string(group, directory_name, "periodic-x", "yes");
	#else
	status = H5LTset_attribute_string(group, directory_name, "periodic-x", "no");
	#endif
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	#ifdef PERIODIC_Y
	status = H5LTset_attribute_string(group, directory_name, "periodic-y", "yes");
	#else
	status = H5LTset_attribute_string(group, directory_name, "periodic-y", "no");
	#endif
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double colloid_diameter=COLLOID_DIAMETER;
	status = H5LTset_attribute_double(group, directory_name, "colloid-diameter", &colloid_diameter, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double colloid_patch_diameter=COLLOID_PATCH_DIAMETER;
	status = H5LTset_attribute_double(group, directory_name, "colloid-patch-diameter", &colloid_patch_diameter, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double colloid_min_bonding_dist=COLLOID_MIN_BONDING_DISTANCE;
	status = H5LTset_attribute_double(group, directory_name, "colloid-min-bonding-distance", &colloid_min_bonding_dist, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double energy_well_depth=ENERGY_WELL_DEPTH;
	status = H5LTset_attribute_double(group, directory_name, "energy-well-depth", &energy_well_depth, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double energy_bond=ENERGY_BOND;
	status = H5LTset_attribute_double(group, directory_name, "energy-bond", &energy_bond, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double temperature=TEMPERATURE;
	status = H5LTset_attribute_double(group, directory_name, "temperature", &temperature, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const unsigned int monte_carlo_steps_main=MONTE_CARLO_STEPS_MAIN;
	status = H5LTset_attribute_uint(group, directory_name, "monte-carlo-steps-main", &monte_carlo_steps_main, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	#ifdef SUBSTRATE_TRIGONAL
	const unsigned int substrate_wells_x=SUBSTRATE_WELLS_X;
	status = H5LTset_attribute_uint(group, directory_name, "substrate-wells-x", &substrate_wells_x, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const unsigned int substrate_wells_y=SUBSTRATE_WELLS_Y;
	status = H5LTset_attribute_uint(group, directory_name, "substrate-wells-y", &substrate_wells_y, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	const double offset_odd=SUBSTRATE_OFFSET_ODD;
	status = H5LTset_attribute_double(group, directory_name, "substrate-offet-odd", &offset_odd, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	#endif
	#ifdef SUBSTRATE_RANDOM
	const unsigned int substrate_number_of_patches=SUBSTRATE_NUMBER_OF_PATCHES;
	status = H5LTset_attribute_uint(group, directory_name, "substrate-number-of-patches", &substrate_number_of_patches, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	#endif
	const double substrate_well_radius=SUBSTRATE_WELL_RADIUS;
	status = H5LTset_attribute_double(group, directory_name, "substrate-well-radius", &substrate_well_radius, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }

	//create simulation frame table
	simframe_size = sizeof( SimulationFrame );
	simframe_offsets[0] = HOFFSET( SimulationFrame, position );
	simframe_offsets[1] = HOFFSET( SimulationFrame, frame_index );
	simframe_offsets[2] = HOFFSET( SimulationFrame, internal_energy );
	simframe_offsets[3] = HOFFSET( SimulationFrame, external_energy );
	simframe_offsets[4] = HOFFSET( SimulationFrame, total_energy );
	simframe_offsets[5] = HOFFSET( SimulationFrame, realtime_seconds );
	simframe_offsets[6] = HOFFSET( SimulationFrame, largest_cluster );
	const char *simframe_field_names[7] = { "position",
	                                        "frame_index",
	                                        "internal_energy",
	                                        "external_energy",
	                                        "total_energy",
						"realtime_seconds",
						"largest_cluster"};
	hid_t simframe_type[7];
	hsize_t flarray_dims[2] = {NUMBER_OF_PARTICLES, 3};
	hid_t flarray_type = H5Tarray_create( H5T_NATIVE_FLOAT, 2, flarray_dims );
	simframe_type[0]=flarray_type;
	simframe_type[1]=H5T_NATIVE_INT;
	simframe_type[2]=H5T_NATIVE_DOUBLE;
	simframe_type[3]=H5T_NATIVE_DOUBLE;
	simframe_type[4]=H5T_NATIVE_DOUBLE;
	simframe_type[5]=H5T_NATIVE_ULONG;
	simframe_type[6]=H5T_NATIVE_INT;

	status = H5TBmake_table("Simulation Frames", //table title
	                        group, //parent node
	                        simulation_frames_location, //path of the dataset
	                        7, //number of fields
	                        0, //number of initial records
	                        simframe_size, //total size of one record
	                        simframe_field_names, //names of the fields
	                        simframe_offsets, //offsets of the fields
	                        simframe_type, //array of types of the fields
	                        20, //chunk size
	                        NULL, //fill data?
	                        5, //compression level (0-9)
	                        NULL); //initial data
	if(status < 0){ printf("> H5Log experienced an error creating the framelog table\n"); return false; }

	//create cluster size table
	cluster_bin_size = sizeof( ClusterSizeBin );
	cluster_bin_offsets[0] = HOFFSET( ClusterSizeBin, frequency );
	cluster_bin_offsets[1] = HOFFSET( ClusterSizeBin, relative_frequency );
	cluster_bin_offsets[2] = HOFFSET( ClusterSizeBin, value );
	const char *cluster_bin_field_names[3] = { "frequency",
	                                        "relative_frequency",
						"value"};
	hid_t cluster_bin_type[3];
	cluster_bin_type[0]=H5T_NATIVE_INT;
	cluster_bin_type[1]=H5T_NATIVE_DOUBLE;
	cluster_bin_type[2]=H5T_NATIVE_INT;

	status = H5TBmake_table("Cluster size distribution", //table title
	                        group, //parent node
	                        cluster_size_distribution_location, //path of the dataset
	                        3, //number of fields
	                        0, //number of initial records
	                        cluster_bin_size, //total size of one record
	                        cluster_bin_field_names, //names of the fields
	                        cluster_bin_offsets, //offsets of the fields
	                        cluster_bin_type, //array of types of the fields
	                        10, //chunk size
	                        NULL, //fill data?
	                        5, //compression level (0-9)
	                        NULL); //initial data
	if(status < 0){ printf("> H5Log experienced an error creating the cluster size distribution table\n"); return false; }
	return true;
}

/**
 * Log the overall acceptance probability, execution time and a few others
 * @param acceptance_probability Average acceptance probability between 1.0 and 0.0.
 * @param execution_time Real execution time in seconds
 * @param max_displacement The maximum displacement used in the simulation. Important if simulation is to be continued.
 * @param max_rotation The maximum rotation used in the simulation. Important if simulation is to be continued.
 * @return Could it be successfully written to the log?
 */
bool h5log_log_statistics(double acceptance_probability, unsigned long execution_time, double max_displacement, double max_rotation){
	printf("> Overall acceptance probability: %3.0f%%\n",acceptance_probability*100);
	herr_t status = H5LTset_attribute_double(group, directory_name, "acceptance-probability", &acceptance_probability, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	printf("> Total execution time: %ld seconds\n",execution_time);
	status = H5LTset_attribute_ulong(group, directory_name, "total-execution-time", &execution_time, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }

	status = H5LTset_attribute_double(group, directory_name, "max-displacement", &max_displacement, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	status = H5LTset_attribute_double(group, directory_name, "max-rotation", &max_rotation, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); return false; }
	return true;
}

/**
 * Log the cluster size distribution
 * @return Was the size successfully logged?
 */
bool h5log_log_cluster_size(void){
	int number_of_bins=0;
	int total_clusters=0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		if( cluster_sizes[i] > 0 ){
			number_of_bins++;
			total_clusters+=cluster_sizes[i];
		}
	}
	ClusterSizeBin *csb=calloc(number_of_bins,sizeof(ClusterSizeBin));
	int j=0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		if( cluster_sizes[i] > 0 ){
			csb[j].frequency=cluster_sizes[i];
			csb[j].relative_frequency=1.0*cluster_sizes[i]/total_clusters;
			csb[j].value=i;
			++j;	
		}
	}
	size_t cluster_bin_sizes[3];
	cluster_bin_sizes[0] = sizeof(csb[0].frequency);
	cluster_bin_sizes[1] = sizeof(csb[0].relative_frequency);
	cluster_bin_sizes[2] = sizeof(csb[0].value);

	herr_t status = H5TBappend_records(group, cluster_size_distribution_location, number_of_bins, cluster_bin_size, cluster_bin_offsets, cluster_bin_sizes, csb);
	if(status < 0){ printf("> H5Log experienced an error loggin a frame\n"); return false; }
	return true;
}


/**
 * Log one frame of the simulation
 * @param particles An array of all particles, frozen. Only energy and position are guaranteed to be correct.
 * @param mc_time The current monte carlo step
 * @param execution_time Real time since start of the simulation in seconds
 * @param largest_cluster The size of the largest cluster in this frame
 * @return Was the frame successfully logged?
 */
bool h5log_log_frame(Colloid *particles, int mc_time, unsigned long execution_time, int largest_cluster){
	SimulationFrame sf[1];
	sf[0].frame_index=mc_time;
	sf[0].realtime_seconds=execution_time;
	sf[0].largest_cluster=largest_cluster;
	sf[0].internal_energy=0.0;
	sf[0].external_energy=0.0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		sf[0].internal_energy+=ENERGY_BOND*particles[i].internal_energy/2.0;
		sf[0].external_energy+=particles[i].external_energy;
		sf[0].position[i][0]=particles[i].position.c.x;
		sf[0].position[i][1]=particles[i].position.c.y;
		sf[0].position[i][2]=particles[i].phi;
	}
	sf[0].total_energy=sf[0].internal_energy+sf[0].external_energy;
	size_t simframe_sizes[7] = { sizeof(sf[0].position),
	                             sizeof(sf[0].frame_index),
	                             sizeof(sf[0].internal_energy),
	                             sizeof(sf[0].external_energy),
	                             sizeof(sf[0].total_energy),
				     sizeof(sf[0].realtime_seconds),
				     sizeof(sf[0].largest_cluster)};

	herr_t status = H5TBappend_records(group, simulation_frames_location, 1, simframe_size, simframe_offsets, simframe_sizes, &sf);
	if(status < 0){ printf("> H5Log experienced an error loggin a frame\n"); return false; }
	return true;
}

/**
 * Close the logfile
 * @return Could the log be successfully written?
 */
bool h5log_close(void){
	herr_t status;
	
	status = H5Gclose(group);
	if(status < 0){ printf("> H5Log experienced an error closing the group\n"); return false; }
	status = H5Fclose(logfile);
	if(status < 0){ printf("> H5Log experienced an error closing the log file\n"); return false; }
	return true;
}

/**
 * Write the patch locations to the logfile
 * 
 * @param centres The locations of the patches
 * @return Was it successfull?
 */
bool h5log_log_substrate(vector2d *centres){
	double buffer[SUBSTRATE_NUMBER_OF_PATCHES*2];
	for(int i=0;i<SUBSTRATE_NUMBER_OF_PATCHES;++i){
		buffer[2*i]=centres[i].c.x;
		buffer[2*i+1]=centres[i].c.y;
	}
	hsize_t dims[2]={SUBSTRATE_NUMBER_OF_PATCHES,2};
	herr_t status = H5LTmake_dataset_double(group,patch_storage_location,2,dims,buffer);
	return status >= 0;
}

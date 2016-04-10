#include <stdio.h>
#include <hdf5.h>
#include <hdf5_hl.h>
#include <x86intrin.h>
#include <stdbool.h>
#include "config.h"
#include "colloid.h"
#include "simulation_frame.h"

hid_t logfile;
hid_t group;
char directory_name[102];
char simulation_frames_location[120];
size_t simframe_size;
size_t simframe_offsets[5];

void h5log_init(void){
	herr_t status;

	logfile = H5Fopen(LOGFILE, H5F_ACC_RDWR, H5P_DEFAULT);
	if(logfile < 0){ //file does not exist
		logfile = H5Fcreate(LOGFILE, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
		if(logfile < 0){
			printf("> H5Log cannot open a file to log this simulation\n");
		}
	}

	//create group
	snprintf(directory_name, sizeof(directory_name), "/%.100s", SIMULATION_SHORT_NAME);
	snprintf(simulation_frames_location, sizeof(simulation_frames_location), "/%.100s/simulation_frames", SIMULATION_SHORT_NAME);

	group = H5Gcreate(logfile, directory_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

	//create group attributes (stores simulation parameters)
	status = H5LTset_attribute_string(group, directory_name, "simulation-name", SIMULATION_NAME);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double size_x=SIZE_X;
	status = H5LTset_attribute_double(group, directory_name, "size-x", &size_x, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double size_y=SIZE_Y;
	status = H5LTset_attribute_double(group, directory_name, "size-y", &size_y, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	#ifdef PERIODIC_X
	status = H5LTset_attribute_string(group, directory_name, "periodic-x", "yes");
	#else
	status = H5LTset_attribute_string(group, directory_name, "periodic-x", "no");
	#endif
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	#ifdef PERIODIC_Y
	status = H5LTset_attribute_string(group, directory_name, "periodic-y", "yes");
	#else
	status = H5LTset_attribute_string(group, directory_name, "periodic-y", "no");
	#endif
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double colloid_diameter=COLLOID_DIAMETER;
	status = H5LTset_attribute_double(group, directory_name, "colloid-diameter", &colloid_diameter, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double colloid_patch_diameter=COLLOID_PATCH_DIAMETER;
	status = H5LTset_attribute_double(group, directory_name, "colloid-patch-diameter", &colloid_patch_diameter, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double colloid_min_bonding_dist=COLLOID_MIN_BONDING_DISTANCE;
	status = H5LTset_attribute_double(group, directory_name, "colloid-min-bonding-distance", &colloid_min_bonding_dist, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double energy_well_depth=ENERGY_WELL_DEPTH;
	status = H5LTset_attribute_double(group, directory_name, "energy-well-depth", &energy_well_depth, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double energy_bond=ENERGY_BOND;
	status = H5LTset_attribute_double(group, directory_name, "energy-bond", &energy_bond, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const double temperature=TEMPERATURE;
	status = H5LTset_attribute_double(group, directory_name, "temperature", &temperature, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }
	const unsigned int monte_carlo_steps_main=MONTE_CARLO_STEPS_MAIN;
	status = H5LTset_attribute_uint(group, directory_name, "monte-carlo-steps-main", &monte_carlo_steps_main, 1);
	if(status < 0){ printf("> H5Log experienced an error setting an attribute\n"); }

	//create simulation frame table
	simframe_size = sizeof( SimulationFrame );
	simframe_offsets[0] = HOFFSET( SimulationFrame, position );
	simframe_offsets[1] = HOFFSET( SimulationFrame, frame_index );
	simframe_offsets[2] = HOFFSET( SimulationFrame, internal_energy );
	simframe_offsets[3] = HOFFSET( SimulationFrame, external_energy );
	simframe_offsets[4] = HOFFSET( SimulationFrame, total_energy );
	const char *simframe_field_names[5] = { "position",
	                                        "frame_index",
	                                        "internal_energy",
	                                        "external_energy",
	                                        "total_energy" };
	hid_t simframe_type[5];
	hsize_t flarray_dims[2] = {NUMBER_OF_PARTICLES, 3};
	hid_t flarray_type = H5Tarray_create( H5T_NATIVE_FLOAT, 2, flarray_dims );
	simframe_type[0]=flarray_type;
	simframe_type[1]=H5T_NATIVE_INT;
	simframe_type[2]=H5T_NATIVE_DOUBLE;
	simframe_type[3]=H5T_NATIVE_DOUBLE;
	simframe_type[4]=H5T_NATIVE_DOUBLE;

	status = H5TBmake_table("Simulation Frames",
	                        group,
	                        simulation_frames_location,
	                        5,
	                        0,
	                        simframe_size,
	                        simframe_field_names,
	                        simframe_offsets,
	                        simframe_type,
	                        10,
	                        NULL,
	                        5,
	                        NULL);
	if(status < 0){ printf("> H5Log experienced an error creating the framelog table\n"); }
}

void h5log_log_frame(Colloid *particles, int mc_time){
	SimulationFrame sf[1];
	sf[0].frame_index=mc_time;
	//TODO: Calculate energy
	sf[0].internal_energy=0.0;
	sf[0].external_energy=0.0;
	sf[0].total_energy=0.0;
	for(int i=0;i<NUMBER_OF_PARTICLES;++i){
		sf[0].position[i][0]=particles[i].position[0];
		sf[0].position[i][1]=particles[i].position[1];
		sf[0].position[i][2]=particles[i].phi;
	}
	size_t simframe_sizes[5] = { sizeof(sf[0].position),
	                             sizeof(sf[0].frame_index),
	                             sizeof(sf[0].internal_energy),
	                             sizeof(sf[0].external_energy),
	                             sizeof(sf[0].total_energy) };

	herr_t status = H5TBappend_records(group, simulation_frames_location, 1, simframe_size, simframe_offsets, simframe_sizes, &sf);
	if(status < 0){ printf("> H5Log experienced an error loggin a frame\n"); }
}

void h5log_close(void){
	herr_t status;
	
	status = H5Gclose(group);
	if(status < 0){ printf("> H5Log experienced an error closing the group\n"); }
	status = H5Fclose(logfile);
	if(status < 0){ printf("> H5Log experienced an error closing the log file\n"); }
}

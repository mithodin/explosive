#define NUMBER_OF_PARTICLES 1

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdbool.h>
#include "simulation_frame.h"

int main(void){
	hid_t file;
	herr_t status;

	file = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	//create simulation frame table
	size_t simframe_size = sizeof( SimulationFrame );
	size_t simframe_offsets[5] = { HOFFSET( SimulationFrame, position ),
								    HOFFSET( SimulationFrame, frame_index ),
								    HOFFSET( SimulationFrame, internal_energy ),
								    HOFFSET( SimulationFrame, external_energy ),
								    HOFFSET( SimulationFrame, total_energy ) };
	const char *simframe_field_names[5] = { "position",
										    "frame_index",
											"internal_energy",
											"external_energy",
											"total_energy" };
	hid_t simframe_type[5];
	hid_t flarray_type = H5Tcopy( H5T_NATIVE_FLOAT );
	H5Tset_size( flarray_type, NUMBER_OF_PARTICLES*3 );
	simframe_type[0]=flarray_type;
	simframe_type[1]=H5T_NATIVE_INT;
	simframe_type[2]=H5T_NATIVE_DOUBLE;
	simframe_type[3]=H5T_NATIVE_DOUBLE;
	simframe_type[4]=H5T_NATIVE_DOUBLE;

	status = H5TBmake_table("Simulation Frames",
						    file,
							"/simframe",
							5,
							100,
							simframe_size,
							simframe_field_names,
							simframe_offsets,
							simframe_type,
							10,
							NULL,
							0,
							NULL);
	
/*	SimulationFrame sf[1] = { { {0.0, 0.0, 0.0}, 0, 10.0, 5.0, 15.0} };
	size_t simframe_sizes[5] = { sizeof(sf[0].position),
								 sizeof(sf[0].frame_index),
								 sizeof(sf[0].internal_energy),
								 sizeof(sf[0].external_energy),
								 sizeof(sf[0].total_energy) };

	H5TBappend_records(file, "/simframe", 1, simframe_size, simframe_offsets, simframe_sizes, &sf); */

	status = H5Fclose(file);
	return 0;
}

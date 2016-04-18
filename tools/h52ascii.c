#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

int main(int argc, char **argv){
	if(argc!=3){
		printf("Usage: h52ascii <filename> <simulation group name>\n");
		return -1;
	}
	char *filename=argv[1];
	char *groupname=argv[2];

	hid_t file,group;

	hsize_t number_of_fields,number_of_records;
	herr_t status;

	//Open file
	if( access( filename, F_OK ) == 0 ){
		file = H5Fopen(filename, H5F_ACC_RDONLY, H5P_DEFAULT);
		if(file < 0){
			printf("Error opening file.\n");
			return -1;
		}
	}else{
		printf("File does not exist.\n");
		return -1;
	}
	
	//Open group
	group = H5Gopen(file,groupname,H5P_DEFAULT);
	if( group < 0 ){
		printf("Error opening group.\n");
		return -1;
	}

	//Read attributes
	unsigned int number_of_particles;
	status = H5LTget_attribute_uint(group,groupname,"number-of-particles",&number_of_particles);
	if( status < 0 ){ printf("Error reading attribute.\n"); return -1; }

	//Get number of fields and records
	status = H5TBget_table_info(group,"simulation_frames",&number_of_fields,&number_of_records);
	if( status < 0 ){
		printf("Error reading table.\n");
		return -1;
	}

	printf("Table has %d fields and %d records\n", (int)number_of_fields, (int)number_of_records);
	
	//Get field info
	char **field_names = calloc(number_of_fields,sizeof(char *));
	for(int i=0;i<number_of_fields;++i){
		field_names[i]=calloc(50,sizeof(char));
	}
	size_t *field_sizes = calloc(number_of_fields,sizeof(size_t));
	size_t *field_offsets = calloc(number_of_fields,sizeof(size_t));
	size_t *type_size = malloc(sizeof(size_t));

	status = H5TBget_field_info(group,"simulation_frames",field_names,field_sizes,field_offsets,type_size);
	if( status < 0){
		printf("Error reading field info.\n");
		return -1;
	}

	size_t position,frame_index,internal_energy,external_energy,total_energy,realtime_seconds,largest_cluster;
	for(int i=0;i<number_of_fields;++i){
		if(strcmp(field_names[i],"position")==0){ position=field_offsets[i]; }
		if(strcmp(field_names[i],"frame_index")==0){ frame_index=field_offsets[i]; }
		if(strcmp(field_names[i],"internal_energy")==0){ internal_energy=field_offsets[i]; }
		if(strcmp(field_names[i],"external_energy")==0){ external_energy=field_offsets[i]; }
		if(strcmp(field_names[i],"total_energy")==0){ total_energy=field_offsets[i]; }
		if(strcmp(field_names[i],"realtime_seconds")==0){ realtime_seconds=field_offsets[i]; }
		if(strcmp(field_names[i],"largest_cluster")==0){ largest_cluster=field_offsets[i]; }
	}

	printf("Size of a frame: %d\n",(int)*type_size);

	//Read frames
	void *table_buffer = calloc(number_of_records,*type_size);
	FILE *file_energy,*file_largest_cluster,*file_movie;
	file_energy=fopen("energy.dat","w");
	if( file_energy == NULL ){
		printf("Error opening energy log file\n");
		return -1;
	}
	file_largest_cluster=fopen("largest_cluster.dat","w");
	if( file_largest_cluster == NULL ){
		printf("Error opening order parameter log file\n");
		return -1;
	}
	file_movie=fopen("movie.xyz","w");
	if( file_movie == NULL ){
		printf("Error opening movie file\n");
		return -1;
	}

	status = H5TBread_table(group, "simulation_frames", *type_size, field_offsets, field_sizes, table_buffer);
	if( status < 0 ){
		printf("Error reading table.\n");
		return -1;
	}

	for(int i=0;i<number_of_records;++i){
		float *positions = table_buffer+(*type_size)*i+position;
		int *frame_index_i = table_buffer+(*type_size)*i+frame_index;
		double *internal_energy_i = table_buffer+(*type_size)*i+internal_energy;
		double *external_energy_i = table_buffer+(*type_size)*i+external_energy;
		double *total_energy_i = table_buffer+(*type_size)*i+total_energy;
		unsigned long *realtime_seconds_i = table_buffer+(*type_size)*i+realtime_seconds;
		int *largest_cluster_i = table_buffer+(*type_size)*i+largest_cluster;
		fprintf(file_largest_cluster,"%d\t%d\n",*frame_index_i,*largest_cluster_i);
		fprintf(file_energy,"%d\t%1.5e\t%1.5e\t%1.5e\n",*frame_index_i,*internal_energy_i,*external_energy_i,*total_energy_i);
		fprintf(file_movie,"%d\n#%d\n",number_of_particles,*frame_index_i);
		for(int p=0;p<number_of_particles;++p){
			fprintf(file_movie,"%3.5f\t%3.5f\t%2.5f\t0\n",positions[3*p],positions[3*p+1],positions[3*p+2]);
		}
	}

	if( fclose(file_movie) != 0 || fclose(file_largest_cluster) != 0 || fclose(file_energy) != 0 ){
		printf("Error closing log files.\n");
		return -1;
	}

	status = H5Gclose(group);
	if( status < 0 ){
		printf("Error closing group\n");
		return -1;
	}

	status = H5Fclose(file);
	if( status < 0 ){
		printf("Error closing file\n");
		return -1;
	}
	return 0;
}

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(int argc, char **argv){
	if(argc!=4){
		printf("Usage: h5clustersize <filename> <simulation group name> <binwidth>\n");
		return -1;
	}
	char *filename=argv[1];
	char *groupname=argv[2];
	int binwidth=atoi(argv[3]);

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

	//Get number of fields and records
	status = H5TBget_table_info(group,"cluster_size",&number_of_fields,&number_of_records);
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
	size_t type_size;

	status = H5TBget_field_info(group,"cluster_size",field_names,field_sizes,field_offsets,&type_size);
	if( status < 0){
		printf("Error reading field info.\n");
		return -1;
	}

	size_t frequency,relative_frequency,value;
	for(int i=0;i<number_of_fields;++i){
		if(strcmp(field_names[i],"frequency")==0){ frequency=field_offsets[i]; }
		if(strcmp(field_names[i],"relative_frequency")==0){ relative_frequency=field_offsets[i]; }
		if(strcmp(field_names[i],"value")==0){ value=field_offsets[i]; }
	}

	printf("Size of a frame: %d\n",(int)type_size);

	//Read frames
	void *table_buffer = calloc(number_of_records,type_size);
	FILE *file_histogram=fopen("cluster_size.dat","w");
	if( file_histogram == NULL ){
		printf("Error opening energy log file\n");
		return -1;
	}

	status = H5TBread_table(group, "cluster_size", type_size, field_offsets, field_sizes, table_buffer);
	if( status < 0 ){
		printf("Error reading table.\n");
		return -1;
	}

	//create histogram
	double *histogram;
	int largest_cluster=0;
	for(int i=0;i<number_of_records;++i){
		int *value_i = table_buffer+type_size*i+value;
		if( *value_i > largest_cluster ){
			largest_cluster=*value_i;
		}
	}
	size_t histogram_num_bins=(size_t)ceil(1.0*largest_cluster/binwidth);
	histogram=calloc(histogram_num_bins,sizeof(double));
	for(int i=0;i<histogram_num_bins;++i){
		histogram[i]=0;
	}
	for(int i=0;i<number_of_records;++i){
		double *relative_frequency_i = table_buffer+type_size*i+relative_frequency;
		int *value_i = table_buffer+type_size*i+value;
		histogram[(int)floor((*value_i-1.0)/binwidth)]+=*relative_frequency_i;
	}
	for(int i=0;i<histogram_num_bins;++i){
		fprintf(file_histogram,"%5.1f\t%3.5f\n",(i+0.5)*binwidth+1.0,histogram[i]);
	}

	if( fclose(file_histogram) != 0 ){
		printf("Error closing log file.\n");
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

	free(field_sizes);
	free(field_offsets);
	free(table_buffer);
	for(int i=0;i<number_of_fields;++i){
		free(field_names[i]);
	}
	free(field_names);
	free(histogram);
	return 0;
}

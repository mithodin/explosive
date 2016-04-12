#include <pthread.h>
#include <semaphore.h>
#include <stdbool.h>
#include <x86intrin.h>
#include <string.h>
#include <math.h>
#include "dSFMT/dSFMT.h"
#include "config.h"
#include "colloid.h"
#include "globals.h"
#include "hdf5_logging.h"

#define LOG_BUFFER_SIZE 5

void *log_logging(void *);

pthread_t log_thread_id;
sem_t log_buffer_full;
sem_t log_buffer_empty;
int log_buffer_write_index;
int log_buffer_read_index;
int *thread_return_status;

Colloid state_buffer[LOG_BUFFER_SIZE][NUMBER_OF_PARTICLES];
int time_index[LOG_BUFFER_SIZE];
int simulation_last_frame;

bool log_init(void){
	if( !h5log_init() ){ return false; }
	sem_init(&log_buffer_full,0,LOG_BUFFER_SIZE);
	sem_init(&log_buffer_empty,0,0);
	log_buffer_write_index=0;
	log_buffer_read_index=0;
	simulation_last_frame=-1;
	thread_return_status=malloc(sizeof(int));
	return pthread_create(&log_thread_id, NULL, log_logging, thread_return_status)==0;
}

void *log_logging(void *arg){
	int *ret_status=(int *)arg;
	*ret_status=-1;
	while(true){
		sem_wait(&log_buffer_empty);
		if( !h5log_log_frame(state_buffer[log_buffer_read_index], time_index[log_buffer_read_index])){
			printf("> error writing to hdf5 log\n");
			*ret_status=-1;
			return (void *)ret_status;
		}
		if(simulation_last_frame==time_index[log_buffer_read_index]){
			*ret_status=0;
			return (void *)ret_status;
		}
		log_buffer_read_index=(log_buffer_read_index+1)%LOG_BUFFER_SIZE;
		sem_post(&log_buffer_full);
	}
	return (void *)ret_status;
}

void log_enqueue(int mc_time, bool simulation_done){
	sem_wait(&log_buffer_full);
	//remember: the pointers in the copy will still point to the active colloids
	memcpy(state_buffer[log_buffer_write_index],particles,sizeof(Colloid)*NUMBER_OF_PARTICLES);
	time_index[log_buffer_write_index]=mc_time;
	if(simulation_done){
		simulation_last_frame=mc_time;
	}
	log_buffer_write_index=(log_buffer_write_index+1)%LOG_BUFFER_SIZE;
	sem_post(&log_buffer_empty);
}

bool log_close(void){
	if( pthread_join(log_thread_id,(void **)&thread_return_status) != 0 || *thread_return_status != 0){
		printf("> Error joining logger thread\n");
		return false;
	}
	return h5log_close();
}

void mkpercent(char *buf, int len, double perc){
	buf[0]='[';
	buf[len-2]=']';
	buf[len-1]='\0';
	for(int i=1;i<len-2;++i){
		double tmp=ceil(perc*(len-3));
		if(i-1<(int)tmp){
			buf[i]='=';
		}else if(i-1==(int)tmp){
			buf[i]='>';
		}else{
			buf[i]=' ';
		}
	}
}

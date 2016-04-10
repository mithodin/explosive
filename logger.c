#include <pthread.h>
#include <semaphore.h>
#include <stdbool.h>
#include <x86intrin.h>
#include <string.h>
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

Colloid state_buffer[LOG_BUFFER_SIZE][NUMBER_OF_PARTICLES];
int time_index[LOG_BUFFER_SIZE];
int simulation_last_frame;

void log_init(void){
	h5log_init();
	sem_init(&log_buffer_full,0,LOG_BUFFER_SIZE);
	sem_init(&log_buffer_empty,0,0);
	log_buffer_write_index=0;
	log_buffer_read_index=0;
	simulation_last_frame=-1;
	pthread_create(&log_thread_id, NULL, log_logging, NULL);
}

void *log_logging(void *arg){
	while(true){
		sem_wait(&log_buffer_empty);
		h5log_log_frame(state_buffer[log_buffer_read_index], time_index[log_buffer_read_index]);
		if(simulation_last_frame==time_index[log_buffer_read_index]){
			return NULL;
		}
		log_buffer_read_index=(log_buffer_read_index+1)%LOG_BUFFER_SIZE;
		sem_post(&log_buffer_full);
	}
	return NULL;
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

void log_close(void){
	pthread_join(log_thread_id,NULL);
	h5log_close();
}

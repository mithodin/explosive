/** @file */
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

pthread_t log_thread_id; /**< Thread ID of the log worker thread */
sem_t log_buffer_full; /**< Semaphore to signal when the log buffer is full */
sem_t log_buffer_empty; /**< Semaphore to signal when the log buffer is empty */
int log_buffer_write_index; /**< Where to write into the log buffer */
int log_buffer_read_index; /**< Where to read from the log buffer */
int *thread_return_status; /**< Signals successfull termination of the log worker thread */

Colloid state_buffer[LOG_BUFFER_SIZE][NUMBER_OF_PARTICLES]; /**< Buffer to store particle information */
int time_index[LOG_BUFFER_SIZE]; /**< Buffer to store mc time index */
unsigned long real_time[LOG_BUFFER_SIZE]; /**< Buffer to store real time */
int simulation_last_frame; /**< Store this so we know when to stop the worker thread */

/**
 * Initialize the logging system including the worker thread
 * @return Was the logging system successfully initialized?
 */
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

/**
 * Log worker thread. We split this so the simulation does not have to wait for I/O operations as much
 * @return *0 on success
 */
void *log_logging(void *arg){
	int *ret_status=(int *)arg;
	*ret_status=-1;
	while(true){
		sem_wait(&log_buffer_empty);
		if( !h5log_log_frame(state_buffer[log_buffer_read_index], time_index[log_buffer_read_index], real_time[log_buffer_read_index])){
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

/**
 * Add a simulation frame to the log buffer. May block when the buffer is full.
 * @param mc_time Current monte carlo step
 * @param simulation_done Is this the last frame?
 * @param runtime Real time since start of the simulation in seconds
 */
void log_enqueue(int mc_time, bool simulation_done, unsigned long runtime){
	sem_wait(&log_buffer_full);
	//remember: the pointers in the copy will still point to the active colloids
	memcpy(state_buffer[log_buffer_write_index],particles,sizeof(Colloid)*NUMBER_OF_PARTICLES);
	time_index[log_buffer_write_index]=mc_time;
	real_time[log_buffer_write_index]=runtime;
	if(simulation_done){
		simulation_last_frame=mc_time;
	}
	log_buffer_write_index=(log_buffer_write_index+1)%LOG_BUFFER_SIZE;
	sem_post(&log_buffer_empty);
}

/**
 * Close the log file and shut down the logging system
 * @return Did the logging system shut down correctly?
 */
bool log_close(void){
	if( pthread_join(log_thread_id,(void **)&thread_return_status) != 0 || *thread_return_status != 0){
		printf("> Error joining logger thread\n");
		return false;
	}
	return h5log_close();
}

/**
 * Log simulation stats (total execution time and average acceptance probability)
 *
 * @param execution_time Total execution time in seconds
 * @param acceptance_probability Average acceptance probability
 * @return Could the stats be successfully written?
 */
bool log_simulation_stats(unsigned long execution_time, double acceptance_probability){
	return h5log_log_final_time(execution_time) && h5log_log_acceptance_probability(acceptance_probability);
}

/**
 * Create a progress bar [===>  ]
 * 
 * @param buf A char buffer of length len to hold the progress bar
 * @param len The length of the buffer. Include one for a null terminator
 * @param perc How much of the progress bar to fill
 */
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

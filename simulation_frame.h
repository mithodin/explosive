/**
 * @file simulation_frame.h
 * @brief Simulation frame log format
 */

/**
 * @brief Data structure to store one simulation frame including some stats
 */
typedef struct _simulation_frame {
	float position[NUMBER_OF_PARTICLES][3];
	int frame_index;
	double internal_energy;
	double external_energy;
	double total_energy;
	unsigned long realtime_seconds;
	int largest_cluster;
} SimulationFrame;

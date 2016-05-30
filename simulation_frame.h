/**
 * @file simulation_frame.h
 * @brief Simulation frame log format
 */

/**
 * @brief Data structure to store one simulation frame including some stats
 */
typedef struct _simulation_frame {
	double position[NUMBER_OF_PARTICLES][3];
	int frame_index;
	double internal_energy;
	double external_energy;
	double total_energy;
	unsigned long realtime_seconds;
	int largest_cluster;
	double bonding_probability;
} SimulationFrame;

/**
 * @brief Data structure to store a cluster size bin
 */
typedef struct _cluster_bin {
	int frequency;
	double relative_frequency;
	int value;
} ClusterSizeBin;

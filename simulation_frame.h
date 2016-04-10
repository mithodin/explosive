typedef struct _simulation_frame {
	float position[NUMBER_OF_PARTICLES][3];
	int frame_index;
	double internal_energy;
	double external_energy;
	double total_energy;
} SimulationFrame;

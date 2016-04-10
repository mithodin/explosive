typedef struct _simulation_frame {
	float position[NUMBER_OF_PARTICLES][3];
	int frame_index;
	double internal_energy;
	double external_energy;
	double total_energy;
} SimulationFrame;

typedef struct _simulation {
	char name[100];
	double size_x;
	double size_y;
	bool periodic_x;
	bool periodic_y;
	double colloid_diameter;
	double colloid_patch_diameter;
	double colloid_min_bonding_distance;
	double energy_well_depth;
	double energy_bond;
	double number_of_particles;
} Simulation;

//give a short name for your simulation. No spaces!
#define SIMULATION_SHORT_NAME "test"
//give a short description of the simulation
#define SIMULATION_NAME "This is a test"
//width and height of the simulation box
#define SIZE_X 2.0
#define SIZE_Y 2.0
//should the box have periodic boundary conditions in x and y? #undef if not.
#define PERIODIC_X
#define PERIODIC_Y
//define the geometry of the colloids
#define COLLOID_DIAMETER 1.0
#define COLLOID_PATCH_DIAMETER 0.11965683746373795115
#define COLLOID_MIN_BONDING_DISTANCE 1.11965683746373795115
//define the energy parameters
#define ENERGY_WELL_DEPTH 0.5
#define ENERGY_BOND 1.0
//temperature in kbT
#define TEMPERATURE 1.0
//how many particles?
#define NUMBER_OF_PARTICLES 2
//how many monte carlo steps?
#define MONTE_CARLO_STEPS_MAIN 10000
//in what interval (in monte carlo steps) should frames be saved to file during the main simulation phase? First and last frame are always saved.
#define LOGGING_INTERVAL 20
#define LOGFILE "simulation_data.h5"
//use these so we don't calculate them millions of times over
#define ONE_THIRD_PI 1.04719755119659763132
#define TWO_THIRDS_PI 2.09439510239319526264

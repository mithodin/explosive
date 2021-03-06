/**
 * @file config.h
 * @brief Configure your simulation
 */

/** give a short name for your simulation. No spaces! */
#define SIMULATION_SHORT_NAME "dev"
/** give a short description of the simulation */
#define SIMULATION_NAME "Development"
/** width of the simulation box */
#define SIZE_X 60.0
/** height of the simulation box */
#define SIZE_Y 60.0
/** should the box have periodic boundary conditions in x? Comment out if not. */
#define PERIODIC_X
/** should the box have periodic boundary conditions in y? Comment out if not. */
#define PERIODIC_Y
/** define the diameter of the colloids */
#define COLLOID_DIAMETER 1.0
/** define the diameter of a patch on the colloid. Do not change. */
#define COLLOID_PATCH_DIAMETER 0.11965683746373795115
/** This should be the sum of COLLOID_DIAMETER and COLLOID_PATCH_DIAMETER */
#define COLLOID_MIN_BONDING_DISTANCE 1.11965683746373795115
/** define the energy of a particle in a well on the substrate */
#define ENERGY_WELL_DEPTH -4.0
/** define the energy of a single bond */
#define ENERGY_BOND 1.0
/** temperature in kbT */
#define TEMPERATURE 0.15
/** how many particles? */
#define NUMBER_OF_PARTICLES 3300
/** how many monte carlo steps? */
#define MONTE_CARLO_STEPS_MAIN 1000000
/** in what interval (in monte carlo steps) should frames be saved to file during the main simulation phase? First and last frame are always saved. */
#define LOGGING_INTERVAL 1000
/** Location of the log file. Can be one global file for all your simulations, the system will handle it! */
#define LOGFILE "./simulation_data.h5"
/** Set substrate pattern. 0 = Random, 1 = Trigonal, 2 = Square */
#define SUBSTRATE_PATTERN 1
/** Define the radius of a well */
#define SUBSTRATE_WELL_RADIUS 8.0
/** Whether or not the potential should be continuous. 0 = Square well, 1 = Continuous 1/r^2, 2 = continuous r^2  */
#define SUBSTRATE_CONTINUOUS 1
/** Set the number of patches */
#define SUBSTRATE_NUMBER_OF_PATCHES 8
/** stepsize for equilibration smoothing*/
#define EQUILIBRATION_SMOOTHING_STEP 500
/** slope (energy change per mc step) cutoff criterion for equilibration */
#define EQUILIBRATION_THRESHOLD_SLOPE 5.0
/** inner loop stepsize for equilibration */
#define EQUILIBRATION_INNER_LOOP 500
/** location for checkpoint after equilibration (maximum 300 chars) */
#define EQUILIBRATION_CHECKPOINT "./"

/** Random initialization of particles? (1: random, 0: square lattice) */
#define PARTICLES_INIT_RANDOM 0

/** define this if an old simulation should be loaded. If the number of particles in this simulation is different from the current one, random particles will be added or removed */
//#define CONTINUE
/** The file to load the old simulation from. If undefined, LOGFILE is used. You HAVE to do this if OLD_LOGFILE is the same as LOGFILE. */
//#define OLD_LOGFILE "./dev.h5"
/** The directory of the old simulation */
//#define OLD_LOGFILE_GROUP "/dev"

#!/usr/bin/env python3

import argparse as ap

parser=ap.ArgumentParser()
parser.add_argument("-n", "--name", help="a short name for the simulation", default="dev")
parser.add_argument("-d" ,"--descr", help="a description of the simulation", default="Development Test")
parser.add_argument("-p", "--particles", help="the number of particles", required=True)
parser.add_argument("-s", "--steps", help="the number of mc steps", type=int, required=True)
parser.add_argument("-l", "--log", help="the logging interval", type=int, default=1000)
arguments=parser.parse_args()

config="""#define SIMULATION_SHORT_NAME \"{sim_name}\"
#define SIMULATION_NAME \"{sim_description}\"
#define SIZE_X 100.0
#define SIZE_Y 86.60254037844386467635
#define PERIODIC_X
#define PERIODIC_Y
#define COLLOID_DIAMETER 1.0
#define COLLOID_PATCH_DIAMETER 0.11965683746373795115
#define COLLOID_MIN_BONDING_DISTANCE 1.11965683746373795115
#define ENERGY_WELL_DEPTH 5.0
#define ENERGY_BOND 1.0
#define TEMPERATURE 0.15
#define NUMBER_OF_PARTICLES {number_of_particles:.0f}
#define MONTE_CARLO_STEPS_MAIN {mc_steps:.0f}
#define LOGGING_INTERVAL {logging_interval:.0f}
#define LOGFILE \"simulation_data.h5\"
#define SUBSTRATE_WELLS_X 5
#define SUBSTRATE_WELLS_Y 5
#define SUBSTRATE_OFFSET_ODD 10
#define SUBSTRATE_WELL_RADIUS 8.0"""

particles_string=arguments.particles.split(",")
particles=[int(i) for i in particles_string]
for particle_number in particles:
    file=open("{p:d}".format(p=particle_number),"w")
    file.write(config.format(sim_name=arguments.name, sim_description=arguments.descr, number_of_particles=particle_number, mc_steps=arguments.steps, logging_interval=arguments.log))
    file.close()

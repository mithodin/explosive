#!/usr/bin/env python3

import argparse as ap
from numpy import floor,log

parser=ap.ArgumentParser()
parser.add_argument("-n", "--name", help="a short name for the simulation", default="dev")
parser.add_argument("-d" ,"--descr", help="a description of the simulation", default="Development Test")
parser.add_argument("-pl", "--particles_lower", help="the lower number of particles", type=int, required=True)
parser.add_argument("-pu", "--particles_upper", help="the upper number of particles", type=int, required=True)
parser.add_argument("-ps", "--particles_step", help="the increment in the number of particles", type=int, default=50)
parser.add_argument("-s", "--steps", help="the number of mc steps", type=int, required=True)
parser.add_argument("-is", "--init_steps", help="the number of mc steps for the first simulation", type=int, default=0)
parser.add_argument("-swx", "--substrate_wells_x", help="the number of wells in x direction", type=int, default=5)
parser.add_argument("-swy", "--substrate_wells_y", help="the number of wells in y direction", type=int, default=6)
parser.add_argument("-l", "--log", help="the logging interval", type=int, default=1000)
arguments=parser.parse_args()

config="""#define SIMULATION_SHORT_NAME \"{sim_name}\"
#define SIMULATION_NAME \"{sim_description}\"
#define SIZE_X 100.0
#define SIZE_Y 103.92304845413263761162
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
#define LOGFILE \"../simulation_data.h5\"
#define SUBSTRATE_WELLS_X {swx:d}
#define SUBSTRATE_WELLS_Y {swy:d}
#define SUBSTRATE_OFFSET_ODD 10
#define SUBSTRATE_WELL_RADIUS 8.0"""

config_continue="""#define SIMULATION_SHORT_NAME \"{sim_name}\"
#define SIMULATION_NAME \"{sim_description}\"
#define SIZE_X 100.0
#define SIZE_Y 103.92304845413263761162
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
#define LOGFILE \"../simulation_data.h5\"
#define SUBSTRATE_WELLS_X {swx:d}
#define SUBSTRATE_WELLS_Y {swy:d}
#define SUBSTRATE_OFFSET_ODD 10
#define SUBSTRATE_WELL_RADIUS 8.0
#define CONTINUE
#define OLD_LOGFILE_GROUP \"{old_sim_name}\""""

num_simulations=int((arguments.particles_upper-arguments.particles_lower)/arguments.particles_step)+1
filename="hy_{{p:0{ns:d}d}}".format(ns=int(floor(log(num_simulations)/log(10))+1))
index_string="-{{n:0{ns:d}d}}".format(ns=int(floor(log(num_simulations)/log(10))+1))
index=0
particle_number=arguments.particles_lower
file=open(filename.format(p=index),"w")
file.write(config.format(sim_name=arguments.name+index_string.format(n=index), sim_description=arguments.descr, number_of_particles=particle_number, mc_steps=arguments.init_steps if arguments.init_steps > 0 else arguments.steps, logging_interval=arguments.log, swx=arguments.substrate_wells_x, swy=arguments.substrate_wells_y))
file.close()
while particle_number <= arguments.particles_upper:
    index+=1
    file=open(filename.format(p=index),"w")
    file.write(config_continue.format(sim_name=arguments.name+index_string.format(n=index), sim_description=arguments.descr, number_of_particles=particle_number, mc_steps=arguments.steps, logging_interval=arguments.log, old_sim_name="/"+arguments.name+index_string.format(n=index-1), swx=arguments.substrate_wells_x, swy=arguments.substrate_wells_y))
    file.close()
    particle_number+=arguments.particles_step

particle_number-=2*arguments.particles_step
while particle_number >= arguments.particles_lower:
    index+=1
    file=open(filename.format(p=index),"w")
    file.write(config_continue.format(sim_name=arguments.name+index_string.format(n=index), sim_description=arguments.descr, number_of_particles=particle_number, mc_steps=arguments.steps, logging_interval=arguments.log, old_sim_name="/"+arguments.name+index_string.format(n=index-1), swx=arguments.substrate_wells_x, swy=arguments.substrate_wells_y))
    file.close()
    particle_number-=arguments.particles_step

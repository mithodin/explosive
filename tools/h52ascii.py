#!/usr/bin/env python3

import tables as tb
import numpy as np
from sys import argv

logfile=tb.open_file(argv[1], 'r')
group=getattr(logfile.root, argv[2])
frames=group.simulation_frames

np.savetxt("internal_energy.dat",frames.cols.internal_energy)
np.savetxt("external_energy.dat",frames.cols.external_energy)
np.savetxt("total_energy.dat",frames.cols.total_energy)

num_particles=group._v_attrs['number-of-particles'][0]
mc_steps=group._v_attrs['monte-carlo-steps-main'][0]
frame_description_format="#frame {{i:0{len:d}d}}\n".format(len=(int)(np.floor(np.log(mc_steps)/np.log(10))+1))

movie=open('movie.xyz','w')
for frame in frames.iterrows():
    movie.write("{np:d}\n".format(np=num_particles))
    movie.write(frame_description_format.format(i=frame[1]))
    for particle in frame[0]:
        movie.write("C {x:1.10f} {y:1.10f} {a:1.10f} 0\n".format(x=particle[0],y=particle[1],a=particle[2]))

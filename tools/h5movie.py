#!/usr/bin/env python3

import tables as tb
import numpy as np
from sys import argv
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, subplot

logfile=tb.open_file(argv[1], 'r')
group=getattr(logfile.root, argv[2])
frames=group.simulation_frames

num_particles=group._v_attrs['number-of-particles'][0]
mc_steps=group._v_attrs['monte-carlo-steps-main'][0]
size_x=group._v_attrs['size-x'][0]
size_y=group._v_attrs['size-y'][0]
col_radius=group._v_attrs['colloid-diameter'][0]/2.0
patch_radius=group._v_attrs['colloid-patch-diameter'][0]/2.0
frame_filename_format="{name}{{i:0{len:d}d}}.png".format(len=(int)(np.floor(np.log(mc_steps)/np.log(10))+1),name=argv[2])

for frame in frames.iterrows():
    fig=plt.figure(1)
    plt.axis([0,size_x,0,size_y])
    ax=fig.add_subplot(1,1,1)
    for particle in frame[0]:
        body=plt.Circle((particle[0],particle[1]), radius=col_radius, color='r', fill=False)
        p0=plt.Circle((particle[0]+0.5*np.cos(particle[2]),particle[1]+0.5*np.sin(particle[2])),radius=patch_radius,color='r', fill=True)
        p1=plt.Circle((particle[0]+0.5*np.cos(2.0/3.0*np.pi+particle[2]),particle[1]+0.5*np.sin(2.0/3.0*np.pi+particle[2])),radius=patch_radius,color='r', fill=True)
        p2=plt.Circle((particle[0]+0.5*np.cos(4.0/3.0*np.pi+particle[2]),particle[1]+0.5*np.sin(4.0/3.0*np.pi+particle[2])),radius=patch_radius,color='r', fill=True)
        ax.add_patch(body)
        ax.add_patch(p0)
        ax.add_patch(p1)
        ax.add_patch(p2)
    plt.savefig(frame_filename_format.format(i=frame[1]))
    plt.clf()

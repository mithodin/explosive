#!/usr/bin/env python3
import tables as tb
import numpy as np
import progressbar as pb
import os
import shutil
from sys import argv
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, subplot
from subprocess import call

def plot_periodic(x,y,a,width,height,px=True,py=True):
    plot=[]
    for x_offset in range(-1 if px else 0,2 if px else 1):
        for y_offset in range(-1 if py else 0,2 if py else 1):
            plot.append(plt.Circle((x+width*x_offset,y+height*y_offset), radius=col_radius, color='r', fill=False))
            plot.append(plt.Circle((x+width*x_offset+0.5*np.cos(a),y+height*y_offset+0.5*np.sin(a)),radius=patch_radius,color='r', fill=True))
            plot.append(plt.Circle((x+width*x_offset+0.5*np.cos(2.0/3.0*np.pi+a),y+height*y_offset+0.5*np.sin(2.0/3.0*np.pi+a)),radius=patch_radius,color='r', fill=True))
            plot.append(plt.Circle((x+width*x_offset+0.5*np.cos(4.0/3.0*np.pi+a),y+height*y_offset+0.5*np.sin(4.0/3.0*np.pi+a)),radius=patch_radius,color='r', fill=True))
    return plot

logfile=tb.open_file(argv[1], 'r')
group=getattr(logfile.root, argv[2])
frames=group.simulation_frames

num_particles=group._v_attrs['number-of-particles'][0]
mc_steps=group._v_attrs['monte-carlo-steps-main'][0]
mc_frames=frames.nrows
size_x=group._v_attrs['size-x'][0]
size_y=group._v_attrs['size-y'][0]
col_radius=group._v_attrs['colloid-diameter'][0]/2.0
patch_radius=group._v_attrs['colloid-patch-diameter'][0]/2.0
frame_filename_format="./tmp/{name}{{i:0{len:d}d}}.png".format(len=(int)(np.floor(np.log(mc_frames)/np.log(10))+1),name=argv[2])
periodic_x=group._v_attrs['periodic-x']==b'yes'
periodic_y=group._v_attrs['periodic-y']==b'yes'

if not os.path.exists("./tmp"):
    os.makedirs("./tmp")
i=0
bar=pb.ProgressBar(max_value=mc_frames)
for frame in bar(frames.iterrows()):
    fig=plt.figure(1)
    plt.axis([0,size_x,0,size_y])
    plt.axes().set_aspect('equal','datalim')
    ax=fig.add_subplot(1,1,1)
    for particle in frame[0]:
        for patch in plot_periodic(particle[0],particle[1],particle[2],size_x,size_y,periodic_x,periodic_y):
            ax.add_patch(patch)
    plt.savefig(frame_filename_format.format(i=i))
    plt.clf()
    i+=1

call(["avconv", "-framerate", "15", "-f", "image2", "-i", "./tmp/test%0{len:d}d.png".format(len=(int)(np.floor(np.log(mc_frames)/np.log(10))+1),name=argv[2]), "-c:v", "h264", "-crf", "1", "test.mp4"])
shutil.rmtree("./tmp")

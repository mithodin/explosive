#!/usr/bin/env python3.3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import Figure, subplot

delta = 0.11965683746373795115;
i=0

data = np.genfromtxt("test.dat")
for line in data:
    if i%10==0 or i==1992:
        fig=plt.figure(1)
        plt.axis([-2,2,-2,2])
        ax=fig.add_subplot(1,1,1)
        circ=plt.Circle((0,0), radius=0.5, color='g', fill=False)
        ax.add_patch(circ)
        c10=plt.Circle((0.5,0),radius=delta/2.0,color='g', fill=True)
        c11=plt.Circle((0.5*np.cos(2.0/3.0*np.pi),0.5*np.sin(2.0/3.0*np.pi)),radius=delta/2.0,color='g', fill=True)
        c12=plt.Circle((0.5*np.cos(4.0/3.0*np.pi),0.5*np.sin(4.0/3.0*np.pi)),radius=delta/2.0,color='g', fill=True)
        ax.add_patch(c10)
        ax.add_patch(c11)
        ax.add_patch(c12)
        c2=plt.Circle((line[5],line[6]), radius=0.5, color='r', fill=False)
        ax.add_patch(c2)
        c20=plt.Circle((line[5]+0.5*np.cos(line[7]),line[6]+0.5*np.sin(line[7])),radius=delta/2.0,color='r', fill=True)
        c21=plt.Circle((line[5]+0.5*np.cos(2.0/3.0*np.pi+line[7]),line[6]+0.5*np.sin(2.0/3.0*np.pi+line[7])),radius=delta/2.0,color='r', fill=True)
        c22=plt.Circle((line[5]+0.5*np.cos(4.0/3.0*np.pi+line[7]),line[6]+0.5*np.sin(4.0/3.0*np.pi+line[7])),radius=delta/2.0,color='r', fill=True)
        ax.add_patch(c20)
        ax.add_patch(c21)
        ax.add_patch(c22)
        plt.savefig("{i:06d}.png".format(i=i))
        plt.clf()
    i+=1

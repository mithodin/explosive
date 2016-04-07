#!/usr/bin/env python3

import numpy as np


def bonded(c1,c2):
    d=distance(c1[0],c1[1],c2[0],c2[1])
    if d < 1.0 or d > 1.11965683746373795115:
        return False
    for i in range(3):
        for j in range(3):
            if distance(c1[0]+0.5*np.cos(c1[2]+2.0/3.0*np.pi*i),c1[1]+0.5*np.sin(c1[2]+2.0/3.0*np.pi*i),c2[0]+0.5*np.cos(c2[2]+2.0/3.0*np.pi*j),c2[1]+0.5*np.sin(c2[2]+2.0/3.0*np.pi*j))<0.11965683746373795115:
                return True
    return False

def distance(x1,y1,x2,y2):
    dx=x1-x2
    dy=y1-y2
    return np.sqrt(dx*dx+dy*dy)

data=np.genfromtxt("test.dat");
for line in data:
    sites=line[0:2]
    c1=line[2:5]
    c2=line[5:8]
    b=line[8]==1
    bb=bonded(c1,c2)
    if ( b and not bb ) or ( not b and bb ):
        print("{c1[0]}\t{c1[1]}\t{c1[2]}\t{c2[0]}\t{c2[1]}\t{c2[2]}".format(c1=c1,c2=c2))

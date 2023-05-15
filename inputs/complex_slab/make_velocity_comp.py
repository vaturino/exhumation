#!/usr/bin/env python3 
import numpy as np
from tqdm import tqdm
import math as m
from scipy.interpolate import griddata
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=600.e3
ymin=0;ymax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; ynum=600
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 0.e3 #crustal length before bending
radius = 350.e3 # radius of curvature
depth_notch = 200.e3
alpha = 70*np.pi/180
gamma = np.pi/2 - alpha
lith = 100.e3

delta = alpha/2
ism = lith/np.cos(delta)


vyr = 5.e-2
yr = 60*60*24*365
v = vyr/yr


x1 = c_len ; 
y1 = ymax - radius
x2 = c_len + (radius)*np.sin(alpha) 
y2 = ymax - radius + radius*(np.cos(alpha))
ip = y2/np.cos(gamma)
x3 = c_len + radius*np.sin(alpha) + ip*np.sin(gamma)
y3 = 0 

dy = ymax - y2
di = dy/np.cos(delta)
dx = di*np.sin(delta)


thick = (lith)/np.sin(alpha)



# create input y-coordinates (with refined region at shallow depth)
ybound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,ymax-ybound,ynum+1-num_refine)
upper_highres = np.linspace(ymax-ybound,ymax,1+num_refine)
yvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(ynum + 1)
C=np.zeros([No_nodes,4],float)

 
ind=0
for j in tqdm(range(ynum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            y = yvals[j]

            C[ind,0] = x
            C[ind,1] = y


            if y >= ymax - lith and x <= (c_len) :
                C[ind,2]=v
                C[ind,3]=0
            # curved portion of crust (i.e. plate boundary)
            if x >= (c_len) and x <= x1  +(x2 - x1)/(y2 - y1)*(y-y1):
                if ((x-x1)**2 + (y-y1)**2) <= radius**2 and ((x-x1)**2 + (y-y1)**2) >= (radius-lith)**2:
                    # Y = y1 - np.sqrt(radius** - (x-x1)**2)
                    C[ind,2]=v*((y-y1)/np.sqrt((x-x1)**2+(y-y1)**2))
                    C[ind,3]=-v*((x-x1)/np.sqrt((x-x1)**2+(y-y1)**2))
            if x >= x1  +(x2 - x1)/(y2 - y1)*(y-y1):
                if x <= x2  +(x3 - x2)/(y3 - y2)*(y-y2) and x >= x2 - thick  +(x3 - x2)/(y3 - y2)*(y-y2):
                    C[ind,2]=v*np.cos(alpha)
                    C[ind,3]=-v*np.sin(alpha)
            
           

            ind=ind+1
        


# write to file
f= open("text_files/base_vel_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
# f.write("# POINTS: %s\n" % (str((xnum+1)*(ynum+1))))
f.write("# Columns: x y velocityx velocityy\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.6E %.6E\n" % (C[k,0], C[k,1], C[k,2],C[k,3]))
f.close() 

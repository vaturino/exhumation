#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm
import math as m

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=400.e3
ymin=0;zmax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; znum=600
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 50.e3 #crustal length before bending
radius = 150.e3 # radius of curvature
depth_notch = 200.e3
alpha = 70*np.pi/180
gamma = 20*np.pi/180
lith = 100.e3

interv = 50.e3

v = 5.e-2
ism = lith/np.cos(alpha/2)
# h = radius - radius*np.cos(alpha)


x1 = c_len 
y1 = zmax
x2 = c_len + interv
y2 = zmax - interv/np.tan(alpha)
x3 = c_len + 2*interv
y3 = y2 - interv/np.tan(alpha - gamma)



# create input y-coordinates (with refined region at shallow depth)
zbound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,zmax-zbound,znum+1-num_refine)
upper_highres = np.linspace(zmax-zbound,zmax,1+num_refine)
zvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(znum + 1)
C=np.zeros([No_nodes,4],float)

 
ind=0
for j in tqdm(range(znum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            z = zvals[j]

            C[ind,0] = x
            C[ind,1] = z


            if z <= zmax and x <= (c_len) :
                C[ind,2]=v
                C[ind,3]=0
            # curved portion of crust (i.e. plate boundary)
            if x >= (c_len) and x <= c_len+interv and z <= y1 + (y2-y1)*((x-x1)/(x2-x1)):
                    C[ind,2]=v*(np.sin(alpha))
                    C[ind,3]=-v*(np.cos(alpha))
            if x >= c_len+interv and z <= y2 + (y3-y2)*((x-x2)/(x3-x2)):
                    C[ind,2]=v*np.sin(alpha - gamma)
                    C[ind,3]=-v*np.cos(alpha - gamma)
            
           

            ind=ind+1
            
           


# write to file
f= open("text_files/base_vel_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
# f.write("# POINTS: %s\n" % (str((xnum+1)*(znum+1))))
f.write("# Columns: x y velocityx velocityy\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.6f %.6f\n" % (C[k,0], C[k,1], C[k,2],C[k,3]))
f.close() 

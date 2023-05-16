#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=400.e3
ymin=0;zmax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1800; znum=900
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 0.e3 #crustal length before bending
radius = 150.e3 # radius of curvature
depth_notch = 100.e3
alpha = 45

vyr = 5.e-2
yr = 60*60*24*365

v = vyr/yr

x2 = c_len + radius*np.sin(alpha)
y2 = zmax + radius*np.cos(alpha) - radius
x3 = 200.e3
y3 = 0




# create input y-coordinates (with refined region at shallow depth)
zbound = 200.e3    # depth of refinement boundary
num_refine = 900   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,zmax-zbound,znum+1-num_refine)
upper_highres = np.linspace(zmax-zbound,zmax,1+num_refine)
zvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(znum + 1)
C=np.zeros([No_nodes,3],float)
 
ind=0
for j in tqdm(range(znum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            z = zvals[j]

            C[ind,0] = x
            C[ind,1] = z

            # flat portion of crust
            if x > 0 and x <= (c_len):
                C[ind,2]=1
            # curved portion of crust (i.e. plate boundary)
            # elif x > (c_len) and x < (c_len + radius*np.sin(alpha)):
            #     x1 = c_len ; 
            #     z1 = zmax - radius
            #     if ((x-x1)**2 + (z-z1)**2) < radius**2 : 
            #         C[ind,2]=1
            elif x > (c_len) and z <= zmax*((x3 - x)/(x3 - c_len)):
                C[ind,2]=1
           

            ind=ind+1
            
           


# write to file
f= open("text_files/base_vel.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
f.write("# Columns: x y velocity\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

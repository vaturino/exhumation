#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
ymin=0;ymax=400.e3
xmin=0;xmax=4*ymax
# number of points (should be a multiple of xmax and ymax)
znum=1200
xnum=4*znum
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 600.e3 #crustal length before bending
radius = 150.e3 # radius of curvature
depth_notch = 100.e3
dip = 45




# create input y-coordinates (with refined region at shallow depth)
zbound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,ymax-zbound,znum+1-num_refine)
upper_highres = np.linspace(ymax-zbound,ymax,1+num_refine)
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
            if x <= (c_len) and z > (ymax - c_thick):
                C[ind,2]=1
            if x > c_len and x <= c_len +ymax and z < (ymax*(1 - (x- c_len)/(ymax))) and z >= (ymax*(1 - (x- c_len)/(ymax))) - c_thick and z >= ymax - 80.e3:
                C[ind,2]=1

            ind=ind+1
            
           


# write to file
f= open("text_files/base_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
f.write("# Columns: x y composition1\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

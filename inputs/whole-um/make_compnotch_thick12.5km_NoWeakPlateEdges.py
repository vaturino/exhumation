#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
ymin=0;ymax=660.e3
xmin=0;xmax=2*ymax
# number of points (should be a multiple of xmax and ymax)
znum=1200
xnum=2*znum
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 200.e3 #crustal length before bending
radius = 250.e3 # radius of curvature
depth_notch = 200.e3
alpha = 45


x2 = c_len + radius*np.sin(alpha)
y2 = ymax + radius*np.cos(alpha) - radius
x3 = 500.e3
y3 = 0
x4 = x2
y4 = ymax + radius*np.cos(alpha) - radius - c_thick/np.cos(alpha)
x5 = x3 - c_thick
y5 = y3



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
            # curved portion of crust (i.e. plate boundary)
            elif x > (c_len) and x < (c_len + radius*np.sin(alpha)):
                x1 = c_len ; 
                z1 = ymax - radius
                if ((x-x1)**2 + (z-z1)**2) < radius**2 and ((x-x1)**2 + (z-z1)**2) >= (radius-c_thick)**2 and z >= ymax - 100.e3: 
                    C[ind,2]=1
            elif x > (c_len + radius*np.sin(alpha)) and z < y2 + ((y3-y2)/(x3-x2))*(x - x2) and z > y4 + ((y5-y4)/(x5-x4))*(x - x4) and z >= ymax - 100.e3:
                C[ind,2]=1

            ind=ind+1
            
           


# write to file
f= open("text_files/base_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
f.write("# Columns: x y composition1\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

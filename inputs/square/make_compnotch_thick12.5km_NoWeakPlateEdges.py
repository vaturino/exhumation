#!/usr/bin/env python3
import numpy
from tqdm import tqdm

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=400.e3
ymin=0;zmax=400.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; znum=1200
# geometry parameters
c_thick = 5.e3 #crustal thickness
depth_notch = 100.e3



# create input y-coordinates (with refined region at shallow depth)
zbound = 200.e3    # depth of refinement boundary
num_refine = 500   # number of grid points in refined upper layer
lower_lowres = numpy.linspace(0,zmax-zbound,znum+1-num_refine)
upper_highres = numpy.linspace(zmax-zbound,zmax,1+num_refine)
zvals = numpy.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(znum + 1)
C=numpy.zeros([No_nodes,3],float)
 
ind=0
for j in tqdm(range(znum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            z = zvals[j]

            C[ind,0] = x
            C[ind,1] = z

            # flat portion of crust
            if x < (xmax) and z <= (zmax - (zmax/xmax)*x) and z >= (zmax - (zmax/xmax)*x) - c_thick and z >= zmax - depth_notch:
                C[ind,2]=1

            ind=ind+1;
            
           


# write to file
f= open("text_files/base_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
f.write("# Columns: x y composition1\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

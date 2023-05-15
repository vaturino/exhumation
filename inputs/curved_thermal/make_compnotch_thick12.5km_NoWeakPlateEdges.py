#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
ymin=0;ymax=660.e3
xmin=0;xmax=3*ymax
# number of points (should be a multiple of xmax and ymax)
ynum=660
xnum=3*ynum

# geometrical parameters, meters
x_gap = 0.e3;         # distance between plate edges and domain sides
x_SP  = 800.e3;        # subducting plate length
depth_notch  = 200e3;   # initial slab depth
radius_outer = 245e3;   # initial slab radius of curvature
slab_dip = 70.;         # slab dip
thick_cr= 10.e3       # overriding plate thickness
ridge_extent = 0.e3; # length of subducting plate that is a "ridge"


# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(ynum + 1)
C=np.zeros([No_nodes,3],float)
 
ind=0
for j in tqdm(range(ynum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)
            y = ymin + j * ((ymax - ymin)/ynum) 

            C[ind,0] = x
            C[ind,1] = y

            # flat portion of crust
            if x <= (x_gap + x_SP - radius_outer) and y >= ymax - thick_cr:
                C[ind,2]=1
            if x > (x_gap + x_SP - radius_outer) and x < (x_gap + x_SP):
                x1 = x_gap + x_SP - radius_outer; 
                y1 = ymax - radius_outer
                if ((x-x1)**2 + (y-y1)**2) < radius_outer**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - thick_cr)**2 and y > (ymax - depth_notch): 
                    angle=np.arctan((y-y1)/(x-x1))
                    if angle > np.radians(90. - slab_dip):
                        ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
                        C[ind,2]=1
                
            ind=ind+1
            
           


# write to file
f= open("text_files/base_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y composition1\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

#!/usr/bin/python3

import sys
import numpy as np
import scipy, scipy.special
import pandas as pd
import matplotlib.pyplot as plt
from tqdm import tqdm

ofile="text_files/part_base.txt"

# box dimensions (i.e. "extent" in ASPECC input)
ymin=800;ymax=900.e3
xmin=0;xmax=6*ymax


# number of cells in input geometry
ynum=3000
xnum=6*ynum


# geometrical parameters, meters
x_gap = 0.e3;         # distance between plate edges and domain sides
x_SP  = xmax/2 - x_gap;        # subducting plate length
depth_notch  = 200e3;   # initial slab depth
radius_outer = 245e3;   # initial slab radius of curvature
slab_dip = 45.;         # slab dip
cthick = 7.5e3
opcthick = 20.e3
oplthick = 80.e3 - opcthick
cutoff = 100.e3
sed = 2.e3
y_coreSP = 25.e3
core_thick = 15.e3
wz_cutoff = 200.e3


# empty array to store geometry
No_nodes= (xnum + 1) * (ynum + 1)
C= pd.DataFrame(np.zeros([No_nodes,3],float))
C.columns = ['x', 'y', 'compo']

 
ind=0
print("writing text file...")
for j in tqdm(range(ynum + 1)): 

	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 
  
		C.iloc[ind,0] = x
		C.iloc[ind,1] = y

		# flat portion of SP
		if  x <= (x_gap + x_SP - radius_outer) and y >= (ymax -cthick-sed) and y <= (ymax):
			C.iloc[ind,2]=1
		

		# curved portion of SP ("notch")
		if x > (x_gap + x_SP - radius_outer) and x < (x_gap + x_SP):
			x1 = x_gap + x_SP - radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) < (radius_outer)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick - sed)**2 and y > (ymax - cutoff): 
				C.iloc[ind,2]= 0
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C.iloc[ind,2]=1	
			
		ind=ind+1;

C = C[C['compo'] != 0].sort_index(ignore_index=True)
num_different_x = len(C['x'].unique())
num_different_y = len(C)/num_different_x
print("Number of x,y which are different:", num_different_x, num_different_y)
print("Number of points: ", len(C))

plt.scatter(C['x'], C['y'], c=C['compo'])
plt.axis('equal')
plt.show()

 
# write to file in ASPECC format
f= open(ofile,"w+")
f.write("# POINTS: %s %s\n" % (num_different_x, int(num_different_y)))
f.write("# Columns: x y compo\n")
C.to_csv(f, sep=" ", header=False, index=False)
f.close()


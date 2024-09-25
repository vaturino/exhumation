#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm

ofile="text_files/comp_base.txt"

# box dimensions (i.e. "extent" in ASPECC input)
ymin=0;ymax=900.e3
xmin=0;xmax=6*ymax


# number of cells in input geometry
ynum=900
xnum=6*ynum


# geometry parameters
x_gap = 400.e3;    # distance between plate edge and box edge
x_SP  = 3000.e3;    # subducting plate length
depth_notch  = 200.e3; # initial depth of crust
radius_outer = 250e3;       # initial radius of curvature of crust
x_box = 50.e3;     # low yielding box width
z_box = 50.e3;     # low yielding box heigth
y_coreSP = 25.e3; # depth of strong core in SP
core_thick = 20.e3; # core thickness
oplthick = 80.e3
wz_cutoff = 200.e3; # cutoff of weak zone
slab_dip = 45.;         # slab dip
cthick = 10e3
cutoff = 100.e3




# empty array to store geometry
No_nodes= (xnum + 1) * (ynum + 1)
C=np.zeros([No_nodes,7],float)
 
ind=0
print("writing text file...")
for j in range(ynum + 1): 

	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 
  
		C[ind,0] = x
		C[ind,1] = y

		### flat portion of all compositions ###
		# SP crust
		if x >= (x_gap) and x <= (x_gap + x_SP - radius_outer) and y >= (ymax -cthick) and y <= (ymax):
			C[ind,2]=1
		# overiding plate (OP) 
		if x >= (x_gap + x_SP) and x <= (xmax - x_gap) and y >= (ymax - oplthick):
			C[ind,5]= 1
		# core of SP
		if x > (x_gap) and x <= (x_gap + x_SP - radius_outer) and y >= (ymax - y_coreSP - core_thick) and y <= (ymax - y_coreSP):
			C[ind,4]=1
		

		### curved portion of all compositions - trench area ###
		if x > (x_gap + x_SP - radius_outer) and x < (x_gap + x_SP):
			x1 = x_gap + x_SP - radius_outer; 
			y1 = ymax - radius_outer;
			# weak zone
			if ((x-x1)**2 + (y-y1)**2) < (radius_outer)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick)**2 and y > (ymax - wz_cutoff): 
				C[ind,6]= 1
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,6]=0
			#subducting plate crust
			if ((x-x1)**2 + (y-y1)**2) < (radius_outer)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick)**2 and y > (ymax - cutoff): 
				C[ind,2]= 0
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,2]=1
			# strong core
			elif ((x-x1)**2 + (y-y1)**2) < radius_outer**2 and ((x-x1)**2 + (y-y1)**2) < (radius_outer-y_coreSP)**2 and ((x-x1)**2 + (y-y1)**2) > (radius_outer-y_coreSP - core_thick)**2 and y > (ymax - cutoff): 
				C[ind,4]= 0
				angle=np.arctan((y-y1)/(x-x1))
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,4]=1
			# overriding plate
			elif ((x-x1)**2 + (y-y1)**2) >= (radius_outer)**2 and  y >= ymax - oplthick:
				C[ind,5]=1
			
		
        # low yielding boxes
		if x > (x_gap - x_box) and x <= (x_gap) and y > (ymax - z_box):  
			C[ind,3]=1
		if x > (xmax - x_gap) and x <= (xmax - x_gap + x_box) and y > (ymax - z_box):
			C[ind,3]=1
		

		
		
			

		
			
		ind=ind+1;
 
# write to file in ASPECC format
f= open(ofile,"w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y compoc compbox compcore compop compwz\n")
for k in range(0,ind):
	f.write("%.4f %.4f %.2f %.2f %.2f %.2f %.2f\n" % (C[k,0],C[k,1],C[k,2],C[k,3],C[k,4],C[k,5],C[k,6]))
f.close() 

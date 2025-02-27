#!/usr/bin/python3

import sys
import numpy as np
import scipy, scipy.special

ofile="text_files/comp_base.txt"

# box dimensions (i.e. "extent" in ASPECC input)
ymin=0;ymax=900.e3
xmin=0;xmax=6*ymax


# number of cells in input geometry
ynum=2400
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
crustal_cutoff = 50.e3
eclogite_cutoff = 150.e3
sed = 2.e3
y_coreSP = 25.e3
core_thick = 20.e3
wz_cutoff = 200.e3
sediment_cutoff = 100.e3


# empty array to store geometry
No_nodes= (xnum + 1) * (ynum + 1)
C=np.zeros([No_nodes,10],float)
 
ind=0
print("writing text file...")
for j in range(ynum + 1): 

	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = ymin + j * ((ymax - ymin)/ynum) 
  
		C[ind,0] = x
		C[ind,1] = y

		# flat portion of SP
		if x >= (x_gap) and x <= (x_gap + x_SP - radius_outer) and y >= (ymax -cthick-sed) and y <= (ymax -sed):
			C[ind,2]=1
		elif x >= (x_gap + x_SP) and x <= (xmax - x_gap) and y >= (ymax -opcthick):
			C[ind,5]= 1
		# overiding plate (OP) - linear geotherm
		elif x >= (x_gap + x_SP) and x <= (xmax - x_gap) and y >= (ymax -opcthick - oplthick):
			C[ind,3]= 1
		elif x >= (x_gap) and x <= (x_gap + x_SP - radius_outer) and y >= (ymax -sed) and y <= (ymax):
			C[ind,4]=1

		# curved portion of SP ("notch")
		if x > (x_gap + x_SP - radius_outer) and x <= (x_gap + x_SP):
			x1 = x_gap + x_SP - radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) <= (radius_outer)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick - sed)**2 and  y > (ymax - wz_cutoff): 
				C[ind,7]= 1
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,7]=0
			if ((x-x1)**2 + (y-y1)**2) <= (radius_outer - sed)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick - sed)**2 and y > (ymax - eclogite_cutoff): 
				C[ind,2]= 0
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,2]=1
			# if ((x-x1)**2 + (y-y1)**2) < (radius_outer - sed)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - cthick - sed)**2 and y <= (ymax - crustal_cutoff) and y >= (ymax - eclogite_cutoff): 
			# 	C[ind,8]= 0
			# 	angle=np.arctan((y-y1)/(x-x1));
			# 	if angle > np.radians(90. - slab_dip):
			# 		ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
			# 		C[ind,8]=1
			if ((x-x1)**2 + (y-y1)**2) <= (radius_outer)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer - sed)**2 and y > (ymax - sediment_cutoff): 
				C[ind,4]= 0
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,4]=1
			elif ((x-x1)**2 + (y-y1)**2) >= (radius_outer)**2 and  y >= (ymax -opcthick):
				C[ind,5]=1
			elif ((x-x1)**2 + (y-y1)**2) >= (radius_outer)**2 and  y >= ymax - (opcthick + oplthick):
				C[ind,3]=1
			
			
		

		
		# flat portion of strong core 
		if x <= (x_SP - radius_outer) and y >= (ymax - y_coreSP - core_thick) and y <= (ymax - y_coreSP):
			C[ind,6]=1
		# curved portion of strong core 
		elif x >= (x_SP - radius_outer):
			x1 = x_gap + x_SP - radius_outer; 
			y1 = ymax - radius_outer;
			if ((x-x1)**2 + (y-y1)**2) < radius_outer**2 and ((x-x1)**2 + (y-y1)**2) <= (radius_outer-y_coreSP)**2 and ((x-x1)**2 + (y-y1)**2) >= (radius_outer-y_coreSP - core_thick)**2 and y > (ymax - depth_notch): 
				C[ind,6]= 0
				angle=np.arctan((y-y1)/(x-x1));
				if angle > np.radians(90. - slab_dip):
					ynotch = radius_outer - np.sqrt((x-x1)**2 + (y-y1)**2)
					C[ind,6]=1
		
			

		
			
		ind=ind+1;
 
# write to file in ASPECC format
f= open(ofile,"w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y compsp compop compsed compopc compcore compwz compecl plstr\n")
for k in range(0,ind):
	f.write("%.4f %.4f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" % (C[k,0],C[k,1],C[k,2],C[k,3],C[k,4],C[k,5],C[k,6],C[k,7], C[k,8], C[k,9]))
f.close() 

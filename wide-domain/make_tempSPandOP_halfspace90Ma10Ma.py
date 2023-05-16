#!/usr/bin/env python3 
import numpy as np, scipy, scipy.special

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
ymin=0;ymax=400.e3
xmin=0;xmax=3*ymax
# number of points (should be a multiple of xmax and ymax)
ynum=1200
xnum=3*ynum

# geometry and temperature parameters
c_thick = 5.e3 #crustal thickness
c_len = 200.e3 #crustal length before bending
radius = 150.e3 # radius of curvature
depth_notch = 200.e3
# subducting and overriding plate ages
age_ma=90;
age=age_ma*1e6*365*24*60*60;
age_op_ma=10;
age_op=age_op_ma*1e6*365*24*60*60;
#--
k = 1e-6                # thermal diffusivity (m^2/s)
Tmax = 1694.5; Tmin = 273.; # min and max temperatures (degrees C)
alpha = 45

x2 = c_len + radius*np.sin(alpha)
y2 = ymax + radius*np.cos(alpha) - radius
x3 = 500.e3
y3 = 0

# create input y-coordinates (with refined region at shallow depth)
ybound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,ymax-ybound,ynum+1-num_refine)
upper_highres = np.linspace(ymax-ybound,ymax,1+num_refine)
yvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, T)
No_nodes= (xnum + 1) * (ynum + 1)
T=np.zeros([No_nodes,3],float)
 
ind=0
for j in range(ynum + 1): 
	for i in range(xnum + 1):

		x = xmin + i * ((xmax - xmin)/xnum)
		y = yvals[j]

		T[ind,0] = x
		T[ind,1] = y
		T[ind,2] = Tmax

		# The thermal structure of all plates is calculated using the ...
		# ... halfspace cooling solution (and specific plate ages)

		# ridge portion of flat subducting plate (i.e. variable age)
		if x >= 0 and x <= (c_len) and y <=ymax:
			erf_term=(ymax-y)/(2*np.sqrt(k*age))
			T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
		# plate boundary region:
		elif x > (c_len) and x < (c_len + radius*np.sin(alpha)) and y<=ymax:
			x1 = c_len; 
			y1 = ymax - radius;
			# subducting plate age below where the crust is (from other script)
			if ((x-x1)**2 + (y-y1)**2) < radius**2: 
				erf_term=(ymax-y)/(2*np.sqrt(k*age))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			# overriding plate age above whether the crust is
			elif ((x-x1)**2 + (y-y1)**2) >= radius**2: 
				erf_term=(ymax-y)/(2*np.sqrt(k*age_op))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
		# flat overriding plate (i.e. uniform age)
		elif x >= (c_len + radius*np.sin(alpha)) and y<=ymax:
			if y < y2 + ((y3-y2)/(x3-x2))*(x - x2):
				erf_term=(ymax-y)/(2*np.sqrt(k*age))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			else:
				erf_term=(ymax-y)/(2*np.sqrt(k*age_op))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))

		ind=ind+1;
 
# write to file
f= open("text_files/temp_base.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y temperature\n")
for k in range(0,ind):
	f.write("%.6f %.6f %.6f\n" % (T[k,0],T[k,1],T[k,2]))
f.close() 


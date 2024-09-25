#!/usr/bin/env python3 
import numpy as np, scipy, scipy.special

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=5400.e3
zmin=0;zmax=900.e3

# number of points (should be a multiple of xmax and ymax)
xnum=5400; znum=900

# geometry parameters
x_gap = 400.e3;    # distance between plate edge and box edge
x_SP  = 3000.e3;    # subducting plate length 
radius = 250e3;
depth_notch  = 200e3;   # initial slab depth
radius = 245e3;   # initial slab radius of curvature
slab_dip = 45.;         # slab dip
opthick = 80e3;        # overriding plate thickness

# thermal parameters
Tmax = 1694.; Tmin = 273.;     # min/max temperatures, K
Tcutoff = 1694.0;               # if T > Tcutoff, T = Tmax
ma=1e6*365*24*60*60;                    # subucting plate age [Ma]
age=90*ma;   # plate age in secs
ridge_extent = 1000.e3; # length of ridge portion of subducting plate
k = 1e-6   

# create input y-coordinates (with refined region at shallow depth)
zbound = 100.e3    # depth of refinement boundary
num_refine = 50   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,zmax-zbound,znum+1-num_refine)
upper_highres = np.linspace(zmax-zbound,zmax,1+num_refine)
zvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, T)
No_nodes= (xnum + 1) * (znum + 1)
T=np.zeros([No_nodes,3],float)
 
ind=0
for j in range(znum + 1): 
		for i in range(xnum + 1):

			x = xmin + i * ((xmax - xmin)/xnum)
			z = zvals[j]

			T[ind,0] = x
			T[ind,1] = z
			T[ind,2] = Tmax

			# The thermal structure of all plates is calculated using the ...
			# ... halfspace cooling solution (and specific plate ages)

			# ridge portion of flat subducting plate (i.e. variable age)
			if x > (x_gap) and x <= (x_gap + ridge_extent):
				age_ridge = (x - x_gap) * (age/ridge_extent)
				erf_term=(zmax-z)/(2*np.sqrt(k*age_ridge))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			# non-ridge portion of flat subducting plate (i.e. uniform age)
			elif x > (x_gap + ridge_extent) and x <= (x_gap + x_SP - radius):
				erf_term=(zmax-z)/(2*np.sqrt(k*age))
				T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			# flat overriding plate (i.e. uniform age)
			elif x >= (x_gap + x_SP) and x <= (xmax - x_gap):
				T[ind,2] = Tmin - (Tmax - Tmin)*(z - zmax)/(opthick)

				# curved portion of SP ("notch")
			if x > (x_gap + x_SP - radius) and x < (x_gap + x_SP):
				x1 = x_gap + x_SP - radius 
				y1 = zmax - radius
				if ((x-x1)**2 + (z-y1)**2) < radius**2 and z > (zmax - depth_notch): 
					erf_term=(zmax-z)/(2*np.sqrt(k*age))
					T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
					angle=np.arctan((z-y1)/(x-x1));
					if angle > np.radians(90. - slab_dip):
						ynotch = radius - np.sqrt((x-x1)**2 + (z-y1)**2)
						erf_term=(ynotch)/(2*np.sqrt(k*age))
						T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
				elif ((x-x1)**2 + (z-y1)**2) >= radius**2:
					T[ind,2] = Tmin - (Tmax - Tmin)*(z - zmax)/(opthick)

			if T[ind,2] > Tcutoff:
				T[ind,2] = Tmax

			# if x > (x_gap + x_SP - radius) and x < (x_gap + x_SP):
			# 	x1 = x_gap + x_SP - radius; 
			# 	z1 = zmax - radius;
			# 	# subducting plate age below where the crust is (from other script)
			# 	if ((x-x1)**2 + (z-z1)**2) < radius**2: 
			# 		erf_term=(zmax-z)/(2*np.sqrt(k*age))
			# 		T[ind,2]='%.5f'  %   (Tmax - (Tmax - Tmin)*scipy.special.erfc(erf_term))
			# 	# overriding plate age above whether the crust is
			# 	elif ((x-x1)**2 + (z-z1)**2) >= radius**2: 
			# 		T[ind,2] = Tmin - (Tmax - Tmin)*(z - zmax)/(opthick)

			ind=ind+1;
 
# write to file
f= open("text_files/temp_base.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1), str(znum+1)))
f.write("# Columns: x y temperature\n")
for k in range(0,ind):
	f.write("%.2f %.2f %.2f\n" % (T[k,0],T[k,1],T[k,2]))
f.close() 


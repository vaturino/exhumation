#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm
import math as m
from scipy.interpolate import griddata
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=400.e3
ymin=0;zmax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; znum=600
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 50.e3 #crustal length before bending
radius = 150.e3 # radius of curvature
depth_notch = 200.e3
alpha = 40*np.pi/180
gamma = np.pi/2 - alpha
lith = 50.e3

delta = alpha/2
v = 5.e-2
ism = lith/np.cos(delta)




x1 = c_len ; 
y1 = zmax - radius
x2 = c_len + (radius)*np.sin(alpha) 
y2 = zmax - radius + radius*(np.cos(alpha))
ip = y2/np.cos(gamma)
x3 = c_len + radius*np.sin(alpha) + ip*np.sin(gamma)
y3 = 0 

dy = zmax - y2
di = dy/np.cos(delta)
dx = di*np.sin(delta)


thick = (lith)/np.sin(alpha)



# create input y-coordinates (with refined region at shallow depth)
zbound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,zmax-zbound,znum+1-num_refine)
upper_highres = np.linspace(zmax-zbound,zmax,1+num_refine)
zvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, z, C)
No_nodes= (xnum + 1) *(znum + 1)
C=np.zeros([No_nodes,4],float)

 
ind=0
for j in tqdm(range(znum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            z = zvals[j]

            C[ind,0] = x
            C[ind,1] = z


            if z >= zmax - lith and x <= (c_len + ism*np.sin(delta)*((z-zmax)/lith)) :
                C[ind,2]=v
                C[ind,3]=0
            # curved portion of crust (i.e. plate boundary)
            if x >= (c_len + ism*np.sin(delta)*((z-zmax)/lith)) and x <= x2 + dx +  (ism*np.sin(delta)*((z-zmax)/lith)):
                if ((x-x1)**2 + (z-y1)**2) <= radius**2 and ((x-x1 + ism*np.sin(delta))**2 + (z-y1 + lith)**2) >= radius**2:
                    # X = x - (x)*np.sin(delta) 
                    X = x - (c_len + dx + ism*np.sin(delta)*((z-zmax)/lith))
                    beta = np.arccos((X - c_len)/((radius)))
                    C[ind,2]=v*(np.sin(beta - delta))
                    C[ind,3]=-v*(np.cos(beta - delta))
            if x >= x2 + dx +  (ism*np.sin(delta)*((z-zmax)/lith)):
                if x <= x2  +(x3 - x2)/(y3 - y2)*(z-y2) and x >= x2 - thick  +(x3 - x2)/(y3 - y2)*(z-y2):
                    C[ind,2]=v*np.cos(alpha)
                    C[ind,3]=-v*np.sin(alpha)
            
           

            ind=ind+1
            
           
print("calculating divergence")
xx = xmax/xnum
yy = zmax/znum
x_low = np.linspace(xmin,xmax,int((xmax-xmin)/xx))
z_low = np.linspace(ymin,zmax,int((zmax-ymin)/yy))
X_low, Z_low = np.meshgrid(x_low,z_low)
# print(X_low.shape, Z_low.shape)
vx= griddata((C[:,0], C[:,0]), C[:,2], (X_low, Z_low), method='cubic')
# vy= griddata((C[:,0], C[:,0]), C[:,3], (X_low, Z_low), method='cubic')

# div = np.zeros((len(vx[:,0]), len(vx[0,:])))

# for j in range(0, len(vx[0,:])-1):
#     for n in range(0, len(vx[:,0])-1):
#         if n == 0 or n == len(vx[:,0]):
#             if j == 0:
#                 div[n, j] = (vx[n, j+1] - vx[n, j])/(dx) + (vy[n+1, j] - vy[n, j])/(dy)
#             elif j == len(vx[0,:]):
#                 div[n, j] = (vx[n, j] - vx[n, j-1])/(dx) + (vy[n, j] - vy[n-1, j])/(dy)
#             else:
#                 div[n, j] = (vx[n, j+1] - vx[n, j-1])/(2*dx)
#         elif j == 0 or j == len(vx[0,:]):
#             if n == 0:
#                 div[n, j] = (vx[n, j+1] - vx[n, j])/(dx) + (vy[n+1, j] - vy[n, j])/(dy)
#             elif n == len(vx[0,:]):
#                 div[n, j] = (vx[n, j] - vx[n, j-1])/(dx) + (vy[n, j] - vy[n-1, j])/(dy)
#             else:
#                 div[n, j] = (vy[n+1, j] - vy[n-1, j])/(2*dy)
#         else:
#             div[n, j] = (vx[n, j+1] - vx[n, j-1])/(2*dx) + (vy[n+1, j] - vy[n-1, j])/(2*dy)


# plt.contourf(X_low, (zmax - Z_low), div, cmap=cm.get_cmap('RdBu_r'),levels=np.linspace(0,1e-10,201),extend='both')
# plt.savefig("test_plots/divergence.png")
# plt.clf()  


# write to file
f= open("text_files/base_vel_comp.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(znum+1)))
# f.write("# POINTS: %s\n" % (str((xnum+1)*(znum+1))))
f.write("# Columns: x y velocityx velocityy\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.6f %.6f\n" % (C[k,0], C[k,1], C[k,2],C[k,3]))
f.close() 

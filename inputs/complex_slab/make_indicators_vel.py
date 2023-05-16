#!/usr/bin/env python3 
import numpy as np
from tqdm import tqdm
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys


# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=600.e3
ymin=0;ymax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; ynum=600
# geometry parameters
c_thick = 5.e3 #crustal thickness
c_len = 0.e3 #crustal length before bending
radius = 350.e3 # radius of curvature
depth_notch = 200.e3
alpha = 70*np.pi/180
gamma = np.pi/2 - alpha
lith = 100.e3

delta =alpha/2

v = 5.e-2
ism = lith/np.cos(delta)
# h = radius - radius*np.cos(alpha)





x1 = c_len ; 
y1 = ymax - radius
x2 = c_len + (radius)*np.sin(alpha) 
y2 = ymax - radius + radius*(np.cos(alpha))
ip = y2/np.cos(gamma)
x3 = c_len + radius*np.sin(alpha) + ip*np.sin(gamma)
y3 = 0 

dy = ymax - y2
di = dy/np.cos(delta)
dx = di*np.sin(delta)


thick = (lith)/np.sin(alpha)





# create input y-coordinates (with refined region at shallow depth)
ybound = 200.e3    # depth of refinement boundary
num_refine = 600   # number of grid points in refined upper layer
lower_lowres = np.linspace(0,ymax-ybound,ynum+1-num_refine)
upper_highres = np.linspace(ymax-ybound,ymax,1+num_refine)
yvals = np.concatenate((lower_lowres, upper_highres[1:]), axis=0)

# create empty array for input file (structure: x, y, y, C)
No_nodes= (xnum + 1) *(ynum + 1)
C=np.zeros([No_nodes,3],float)
 
ind=0
for j in tqdm(range(ynum + 1)): 
        for i in range(xnum + 1):

            x = xmin + i * ((xmax - xmin)/xnum)

            y = yvals[j]

            C[ind,0] = x
            C[ind,1] = y

            if y >= ymax - lith and x <= (c_len) :
                C[ind,2]=1
            # curved portion of crust (i.e. plate boundary)
            if x >= (c_len) and x <= x1  +(x2 - x1)/(y2 - y1)*(y-y1):
                if ((x-x1)**2 + (y-y1)**2) <= radius**2 and ((x-x1)**2 + (y-y1)**2) >= (radius-lith)**2:
                        C[ind,2]=1
            if x >= x1  +(x2 - x1)/(y2 - y1)*(y-y1):
                if x <= x2  +(x3 - x2)/(y3 - y2)*(y-y2) and x >= x2 - thick  +(x3 - x2)/(y3 - y2)*(y-y2):
                    C[ind,2]=1
            
            

            ind=ind+1

    
           


# write to file
f= open("text_files/base_vel.txt","w+")
f.write("# POINTS: %s %s\n" % (str(xnum+1),str(ynum+1)))
f.write("# Columns: x y velocity\n")
for k in range(0,ind):
    f.write("%.6f %.6f %.2f\n" % (C[k,0],C[k,1],C[k,2]))
f.close() 

fig=plt.figure()
gs=GridSpec(1,1) 

ax1=fig.add_subplot(gs[0,0], aspect=1)
input_plot = ax1.scatter(C[:,0]/1.e3, (ymax-C[:,1])/1.e3, c=C[:,2],cmap='bwr',vmin=min(C[:,2]),vmax=max(C[:,2]),s=0.25,lw=0)
# ax1.scatter(x2/1.e3, (ymax -y2)/1.e3, c ='k', s = 0.5)
ax1.set_ylim([(ymax)/1.e3,0])   
ax1.set_xlim([0,(xmax)/1.e3])
ax1.tick_params(labelsize=6)
cbar = plt.colorbar(input_plot, cax = fig.add_axes([0.7, 0.4, 0.15, 0.015]), orientation='horizontal',ticks=[min(C[:,2]), 0.5*(min(C[:,2])+max(C[:,2])), max(C[:,2])], ticklocation = 'top')
cbar.ax.tick_params(labelsize=6)


plot_name = "test_plots/base_vel.png"
plt.savefig(plot_name, bbox_inches='tight', format='png', dpi=500)
plt.clf()

#!/usr/bin/env python3 

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

file_name='base_vel_comp'	            
cmin=float(sys.argv[1])  # minimum value for color bar
cmax=float(sys.argv[2])  # maximum value for color bar

xmax=400.e3; zmax=200.e3

file_path = ''.join(['text_files/',str(file_name),'.txt'])
plot_name = ''.join(['test_plots/',str(file_name),'.png'])
file 	  = np.loadtxt(file_path,skiprows=2)

fig=plt.figure(figsize = (20,15))
gs=GridSpec(3,1) 

print("input file...")
ax1=fig.add_subplot(gs[0,0], aspect=1)
input_plot1 = ax1.scatter(file[:,0]/1.e3, (zmax-file[:,1])/1.e3, c=np.sqrt(file[:,2]**2 + file[:,3]**2),cmap='bwr',vmin=0,vmax=cmax,s=0.25,lw=0)
ax1.set_ylim([(zmax)/1.e3,0])   
ax1.set_xlim([0,(xmax)/1.e3])
ax1.tick_params(labelsize=6)
cbar1 = plt.colorbar(input_plot1, cax = fig.add_axes([0.7, 0.4, 0.15, 0.015]), orientation='horizontal',ticks=[cmin, 0.5*(cmin+cmax), cmax], ticklocation = 'top')
cbar1.ax.tick_params(labelsize=6)

ax2=fig.add_subplot(gs[1,0], aspect=1)
input_plot2 = ax2.scatter(file[:,0]/1.e3, (zmax-file[:,1])/1.e3, c=file[:,2],cmap='bwr',vmin=cmin,vmax=cmax,s=0.25,lw=0)
ax2.set_ylim([(zmax)/1.e3,0])   
ax2.set_xlim([0,(xmax)/1.e3])
ax2.tick_params(labelsize=6)
# cbar2 = plt.colorbar(input_plot2, cax = fig.add_axes([0.7, 0.4, 0.15, 0.015]), orientation='horizontal',ticks=[cmin, 0.5*(cmin+cmax), cmax], ticklocation = 'top')
# cbar2.ax.tick_params(labelsize=6)

ax3=fig.add_subplot(gs[2,0], aspect=1)
input_plot3 = ax3.scatter(file[:,0]/1.e3, (zmax-file[:,1])/1.e3, c=abs(file[:,3]),cmap='bwr',vmin=cmin,vmax=cmax,s=0.25,lw=0)
ax3.set_ylim([(zmax)/1.e3,0])   
ax3.set_xlim([0,(xmax)/1.e3])
ax3.tick_params(labelsize=6)
# cbar3 = plt.colorbar(input_plot3, cax = fig.add_axes([0.7, 0.4, 0.15, 0.015]), orientation='horizontal',ticks=[cmin, 0.5*(cmin+cmax), cmax], ticklocation = 'top')
# cbar3.ax.tick_params(labelsize=6)


print("saving plot to %s..." % plot_name)
plt.savefig(plot_name, bbox_inches='tight', format='png', dpi=1000)
plt.clf()

   

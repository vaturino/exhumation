#!/usr/bin/python3

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import sys

file_name=str(sys.argv[1])	            
cmin=float(sys.argv[2])  # minimum value for color bar
cmax=float(sys.argv[3])  # maximum value for color bar

zmax=900.e3; xmax=6*zmax

file_path = ''.join(['text_files/',str(file_name),'.txt'])
plot_name = ''.join(['test_plots/',str(file_name),'.eps'])
file 	  = np.loadtxt(file_path,skiprows=2)

fig=plt.figure()
gs=GridSpec(1,1) 

print("input file...")
ax1=fig.add_subplot(gs[0,0], aspect=1)
if file_name == 'comp_base':
    input_plot = ax1.scatter(file[:,0]/1.e3, (zmax-file[:,1])/1.e3, c=file[:,2]+2*file[:,3]+3*file[:,4]+4*file[:,5]+5*file[:,6]+6*file[:,7]+7*file[:,8] + 8*file[:,9],cmap='bwr',vmin=cmin,vmax=cmax,s=0.25,lw=0)
else:
    input_plot = ax1.scatter(file[:,0]/1.e3, (zmax-file[:,1])/1.e3, c = file[:,2],cmap='bwr',vmin=cmin,vmax=cmax,s=0.25,lw=0)
ax1.set_ylim([(zmax)/1.e3,0])   
ax1.set_xlim([0,(xmax)/1.e3])
ax1.tick_params(labelsize=6)
cbar = plt.colorbar(input_plot, cax = fig.add_axes([0.7, 0.4, 0.15, 0.015]), orientation='horizontal',ticks=[cmin, 0.5*(cmin+cmax), cmax], ticklocation = 'top')
cbar.ax.tick_params(labelsize=6)

print("saving plot to %s..." % plot_name)
plt.savefig(plot_name, bbox_inches='tight', format='eps', dpi=1000)
plt.clf()

   

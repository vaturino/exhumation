#!/usr/bin/env python3 
import numpy as np
from tqdm import tqdm
import math as m
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from scipy.interpolate import griddata

xmin=0;xmax=600.e3
ymin=0;ymax=200.e3
xnum=1200; ynum=600

df = pd.read_csv("text_files/base_vel_comp.txt", sep="\s+", skiprows=2, header = None)
df.columns=['x','y','vx','vy']
df.x = df.x
df.y = ymax - df.y

dx = xmax/xnum
dy = ymax/ynum

x_low = np.linspace(xmin,xmax,int((xmax-xmin)/dx))
y_low = np.linspace(ymin,ymax,int((ymax-ymin)/dy))
X_low, Y_low = np.meshgrid(x_low,y_low)

vx= griddata((df.x, df.y), df.vx, (X_low, Y_low), method='cubic')
vy= griddata((df.x, df.y), df.vy, (X_low, Y_low), method='cubic')

div = np.zeros((len(vx[:,0]), len(vx[0,:])))

for j in range(0, len(vx[0,:])-1):
    for n in range(0, len(vx[:,0])-1):
        if n == 0 or n == len(vx[:,0]):
            if j == 0:
                div[n, j] = (vx[n, j+1] - vx[n, j])/(dx) + (vy[n+1, j] - vy[n, j])/(dy)
            elif j == len(vx[0,:]):
                div[n, j] = (vx[n, j] - vx[n, j-1])/(dx) + (vy[n, j] - vy[n-1, j])/(dy)
            else:
                div[n, j] = (vx[n, j+1] - vx[n, j-1])/(2*dx)
        elif j == 0 or j == len(vx[0,:]):
            if n == 0:
                div[n, j] = (vx[n, j+1] - vx[n, j])/(dx) + (vy[n+1, j] - vy[n, j])/(dy)
            elif n == len(vx[:,0]):
                div[n, j] = (vx[n, j] - vx[n, j-1])/(dx) + (vy[n, j] - vy[n-1, j])/(dy)
            else:
                div[n, j] = (vy[n+1, j] - vy[n-1, j])/(2*dy)
        else:
            div[n, j] = (vx[n, j+1] - vx[n, j-1])/(2*dx) + (vy[n+1, j] - vy[n-1, j])/(2*dy)


plt.contourf(X_low/1.e3, (ymax - Y_low)/1.e3, div, cmap=cm.get_cmap('RdBu_r'),levels=np.linspace(0,1e-10,501),extend='both')
ax=plt.gca()
plt.xlim(0, 400)
plt.ylim(0, 200)
ax.set_aspect(1)
plt.colorbar()
plt.savefig("test_plots/divergence.png", dpi=500)
plt.clf()


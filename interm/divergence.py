#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm
import math as m
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec



def divergence(f,sp):
    """ Computes divergence of vector field 
    f: array -> vector field components [Fx,Fy,Fz,...]
    sp: array -> spacing between points in respecitve directions [spx, spy,spz,...]
    """
    num_dims = len(f)
    return np.ufunc.reduce(np.add, [np.gradient(f[i], sp[i], axis=i) for i in range(num_dims)])

a = []
for N in range(20,100):
    # Number of points (NxN)
    # = 20
    # Boundaries
    ymin = -2.; ymax = 2.
    xmin = -2.; xmax = 2.
    
    
    # Divergence function
    def divergence(f,sp):
        num_dims = len(f)
        return np.ufunc.reduce(np.add, [np.gradient(f[i], sp[i], axis=i) for i in range(num_dims)])

    
    # Create Meshgrid
    x = np.linspace(xmin,xmax, N)
    y = np.linspace(ymin,ymax, N)
    xx, yy = np.meshgrid(x, y)
    
    
    # Define 2D Vector Field
    Fx  = np.cos(xx + 2*yy)
    Fy  = np.sin(xx - 2*yy)
    
    F = np.array([Fx, Fy])
    # Compute Divergence
    sp_x = np.diff(x)[0]
    sp_y = np.diff(y)[0]
    sp = [sp_x, sp_y]
    g = divergence(F, sp)
    
    print("Max: ", np.max(g.flatten()))
    a.append(np.max(g.flatten()))
plt.plot(a)
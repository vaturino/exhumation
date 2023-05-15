#!/usr/bin/env python3
import numpy as np
from tqdm import tqdm
import math as m
import pandas as pd
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec


# box dimensions (make sure consistent with "X extent" and "Y extent" in ASPECT input)
xmin=0;xmax=400.e3
ymin=0;zmax=200.e3
# number of points (should be a multiple of xmax and ymax)
xnum=1200; znum=600

df = pd.read_csv("text_files/base_vel_comp.txt", sep="\s+", skiprows=2, header = None)
df.columns=['x','y','vx','vy']
df.x = df.x/1.e3
df.y = (zmax-df.y)/1.e3

dx = xmax/xnum*1e-3
dy = zmax/znum*1e-3
d = [dx, dy]

F = [df.vx, df.vy]

# print(len(F))

grad = np.gradient(np.array(df.vx, df.vy), dx, edge_order=2)
div = np.sum(grad)

# dim = len(F)
# div = np.ufunc.reduce(np.add, [np.gradient(F)])
# for i in range(len(div)):
#     if div[i].any() != 0:
#         print(div[i])


# div = 0
# for i in range(len(grad)):
#     div = div + grad[i]

    
print(div)

# # write to file
# f= open("text_files/gradient.txt","w+")
# for k in range(0,len(grad)):
#     f.write("%.2E \n" % (grad[k]))
# f.close() 


















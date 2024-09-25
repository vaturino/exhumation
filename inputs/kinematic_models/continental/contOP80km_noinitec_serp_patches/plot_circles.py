#!/usr/bin/python3

import matplotlib.pyplot as plt
import numpy as np

# Parameters for Circle C1
radius = 183.125e3
center_C1 = (0, radius - 5.e3)

# Parameters for spacing between centers of circles
# spacing = 15.e3 + 42.650e3 + np.sqrt(radius**2 - center_C1[1]**2)
spacing = 57.5e3 + np.sqrt(radius**2 - center_C1[1]**2)

n = 7  # Number of circles
num_points = 100
circles = np.zeros((2, num_points, n))

# Define function to create circle data
def create_circle(center, radius, num_points=100):
    theta = np.linspace(0, 2 * np.pi, num_points)
    x = center[0] + radius * np.cos(theta)
    y = center[1] + radius * np.sin(theta)
    return x, y

# Create grid of points
grid_resolution = 500
x_vals = np.linspace(-radius, n * spacing + radius, grid_resolution)
y_vals = np.linspace(center_C1[1] - radius, center_C1[1] + radius, grid_resolution)
X, Y = np.meshgrid(x_vals, y_vals)

# Find points inside any of the circles
inside = np.zeros_like(X, dtype=bool)
for i in range(n):
    circles[0, :, i], circles[1, :, i] = create_circle((i * spacing, center_C1[1]), radius, num_points)
    center = (i * spacing, center_C1[1])
    inside |= (X - center[0])**2 + (Y - center[1])**2 <= radius**2



# Plot the circles and the points inside them
plt.plot(circles[0, :, :], circles[1, :, :], label='Circles')
plt.scatter(X[inside], Y[inside], color='red', s=0.1, label='Points Inside')
plt.gca().set_aspect('equal', adjustable='box')
plt.ylim(-10.e3, 0.e3)
plt.xlabel('x')
plt.ylabel('y')
plt.legend()
plt.grid(True)
plt.show()

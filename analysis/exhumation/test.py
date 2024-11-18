#!/usr/bin/env python3
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

# Create a 2x2 grid of subplots
fig, axs = plt.subplots(2, 2, figsize=(8, 8))

# Example data
for i in range(2):
    for j in range(2):
        axs[i, j].plot([0, 1], [0, 1], label=f'Plot ({i+1}, {j+1})')
        axs[i, j].legend()

# Margins for the rectangles
margin = 0.05  # General margin
margin_right = 0.01

# --- Rectangle for the first column ---
# Get positions of subplots [0, 0] and [1, 0]
bbox1 = axs[0, 0].get_position()
bbox2 = axs[1, 0].get_position()

# Calculate the bounds
x0_col1 = bbox1.x0 - margin  # Left edge
y0_col1 = bbox2.y0 - margin  # Bottom edge
width_col1 = bbox1.width + margin + margin_right  # Width
height_col1 = (bbox1.y1 - bbox2.y0) + margin * 2  # Height

# Add rectangle for the first column
rect_col1 = Rectangle(
    (x0_col1, y0_col1), width_col1, height_col1, transform=fig.transFigure,
    color='red', fill=False, lw=2
)
fig.patches.append(rect_col1)

# --- Rectangle for the second column ---
# Get positions of subplots [0, 1] and [1, 1]
bbox3 = axs[0, 1].get_position()
bbox4 = axs[1, 1].get_position()

# Calculate the bounds
x0_col2 = bbox3.x0 - margin  # Left edge
y0_col2 = bbox4.y0 - margin  # Bottom edge
width_col2 = bbox3.width + margin + margin_right  # Width
height_col2 = (bbox3.y1 - bbox4.y0) + margin * 2  # Height

# Add rectangle for the second column
rect_col2 = Rectangle(
    (x0_col2, y0_col2), width_col2, height_col2, transform=fig.transFigure,
    color='blue', fill=False, lw=2
)
fig.patches.append(rect_col2)

# Show the figure
plt.show()

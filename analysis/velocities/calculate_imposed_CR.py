#! /usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#get convergence rate data: 
cr = pd.read_csv("/home/vturino/PhD/projects/exhumation/plots/single_models/ocean_continent_subduction_fields_longSP_visc24/txt_files/2D_v.txt", sep="\s+")
cr["conv_rate"].iloc[0] = np.nan
cr["conv_rate"] = cr["conv_rate"]

# Define the time ranges for each equation
t1 = np.linspace(0, 7.17225e6, 100)
t2 = np.linspace(7.17225e6, 13.7765e6, 100)
t3 = np.linspace(13.7765e6, 17.7713e6, 100)
t4 = np.linspace(17.7713e6, 35e6, 100)  # Changed to 35e6

# Original equations
m3 = -1.11833e-08
b3 = 0.231524
v1_sp = 1.00756e-09 * t1 + 0.020515
v2_sp = 7.52792e-09 * t2 - 0.026251
v3_final = m3 * t3[-1] + b3  # Calculate the last value of v3
v4_value = 0.034222  # Constant value for v4
v4_sp = np.full_like(t4, v4_value)  # v4 is constant

# Adjust v3 to ensure v3[-1] == v4[0]
v3 = m3 * t3 + b3
v3[-1] = v4_sp[0]  # Ensure v3's last value matches v4's first value

# Define an array of percentages for the linear decrease after t4
percentages = [0.0, 0.1, 0.5, 1.0]  # 0%, 10%, 50%, 100%
t5 = 50e6

# Prepare the output file
output_file = "/home/vturino/PhD/projects/exhumation/plots/imposed_velocities/velocity_equations.txt"
with open(output_file, "w") as file:
    # Write the equations for v1 to v4
    equations = [
        (1.00756e-09, 0.020515, 0, 7.17225e6),
        (7.52792e-09, -0.026251, 7.17225e6, 13.7765e6),
        (m3, b3, 13.7765e6, 17.7713e6),
        (0.0, v4_value, 17.7713e6, 35e6)  # v4 is constant at v4
    ]

    for i, (m, b, t_prev, t_cur) in enumerate(equations, start=1):
        equation_str = f"v_{i} = {m:.6e} * t + {b:.6f}, for t in [{t_prev:.2e}, {t_cur:.2e}]\n"
        file.write(equation_str)

    # Write the equations for the linear decrease based on percentages
    for i, percentage in enumerate(percentages, start=1):
        v5 = percentage * v4_value
        m = (v5 - v4_value) / (t5 - 35e6)  # Slope for the linear decrease
        b = v4_value - m * 35e6            # Intercept for the linear decrease
        
        percentage_str = f"v_5 = {m:.6e} * t + {b:.6f}, for t in [35e6, {t5:.2e}], and v_5 = {percentage*100:.0f}% of v_4\n"
        file.write(percentage_str)

# Plot the equations
plt.figure(figsize=(10, 6))

# Plot the first four segments
plt.plot(cr["time"], cr["conv_rate"]/1e2, label="Dynamic Convergence Rate", color='black')
plt.plot(t1, v1_sp, label=r"$v_{1,sp} = 1.00756 \times 10^{-9} t + 0.020515$", color='b')
plt.plot(t2, v2_sp, label=r"$v_{2,sp} = 7.52792 \times 10^{-9} t - 0.026251$", color='g')
plt.plot(t3, v3, label=r"$v_{3,sp} = -1.11833 \times 10^{-8} t + 0.231524$", color='r')
plt.plot(t4, v4_sp, label=r"$v_{4,sp} = 0.034222$", color='orange')

# After v4, plot linear decreases for each percentage
colors = plt.cm.Dark2(np.arange(len(percentages)))

for i, percentage in enumerate(percentages):
    v5 = percentage * v4_value  # Final velocity value for this percentage
    m = (v5 - v4_value) / (t5 - 35e6)  # Slope for the linear decrease
    t_after = np.linspace(35e6, t5, 100)
    v_after = m * (t_after - 35e6) + v4_value
    plt.plot(t_after, v_after, label=f"Linear decrease to {percentage*100:.0f}% of v4", color=colors[i])

# Add titles and labels
plt.title("Plot of v_sp Equations over Time with Linear Decrease for Multiple Percentages")
plt.xlabel("Time (t)")
plt.ylabel("v_sp(t)")
plt.grid(True)
plt.legend()

# Save the plot
plot_file = "/home/vturino/PhD/projects/exhumation/plots/imposed_velocities/"
plt.savefig(f'{plot_file}/velocity_plot.eps', dpi=1000, bbox_inches='tight', format='eps')
plt.close()

print(f"Equations written to {output_file}")
print(f"Plot saved as {plot_file}")



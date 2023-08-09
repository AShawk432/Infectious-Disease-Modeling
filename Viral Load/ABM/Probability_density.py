import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the viral load data from the CSV file
data = pd.read_csv('viral_load.csv', header=None)
# Prepare data for plotting

# Multiply the viral load data by 10 and round it
viral_loads = np.round(data.values * 10)

# Initialize time steps based on number of columns
time_steps = np.arange(viral_loads.shape[1])  # Assume time steps start from 1

# Create a 3D surface plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Define the number of bins for the viral load histogram
num_bins = 10

# Compute the minimum rounded viral load value divided by 10
min_viral_load = np.min(viral_loads) / 10

# Create arrays to store the bin edges and heights for each time step
bin_edges = np.linspace(0, np.max(viral_loads), num_bins + 1)
bin_heights = np.zeros((num_bins, len(time_steps)))

# Iterate over each time step
for i, time_step in enumerate(time_steps):
    # Compute the histogram for the current time step
    hist, _ = np.histogram(viral_loads[:, i], bins=bin_edges)

    # Store the bin heights
    bin_heights[:, i] = hist / len(viral_loads)   # Compute normalized probability density

# Create a meshgrid for the bin edges and time steps
X, Y = np.meshgrid(bin_edges[:-1], time_steps)
Z = bin_heights.T


# Create the surface plot
surf = ax.plot_surface(X, Y, Z, cmap='viridis', edgecolor='none')

# Set the labels and title
ax.set_xlabel('Viral Load')
ax.set_ylabel('Time Steps')
ax.set_zlabel('Probability Density')
ax.set_title('3D Probability Distribution of Viral Load')

# Add a colorbar
# fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.15)
# Show the plot
plt.show()

# Create a 2D plot with colorbar
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
im = ax2.imshow(Z, cmap='viridis', aspect='auto', extent=[0, 1, time_steps[1], time_steps[0]])
fig2.colorbar(im, ax=ax2)

# Set the labels and title for the 2D plot
ax2.set_xlabel('Viral Load')
ax2.set_ylabel('Time Steps')
ax2.set_title('2D Probability Distribution of Viral Load')
# Reverse the time steps for the 2D plot
ax2.invert_yaxis()

# Show the plots
plt.show()


# # Create an empty list to store the mean viral load for each time step
# mean_viral_loads = []
#
# # Compute the mean viral load for each time step
# for i, time_step in enumerate(time_steps):
#     non_zero_values = viral_loads[:, i][viral_loads[:, i] != 0]  # Select non-zero values at the current time step
#     if len(non_zero_values) > 0:  # Check if there are non-zero values
#         mean_viral_load = np.mean(non_zero_values)
#         mean_viral_loads.append(mean_viral_load)
#
# # Plot the viral load over time
# plt.plot(time_steps[:len(mean_viral_loads)], mean_viral_loads, marker='o')  # Use sliced time_steps for non-zero mean values
# plt.xlabel('Time Steps')
# plt.ylabel('Viral Load')
# plt.title('Non-Zero Mean Viral Load Over Time')
# plt.grid(True)
#
# # Show the plot
# plt.show()

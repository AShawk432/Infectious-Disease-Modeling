import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Specify the directory containing the Excel files
data_directory = 'Viral_Load_Data'

# Create a list of file names in the directory
file_list = os.listdir(data_directory)

# Create a directory for save the plots
output_directory = 'VL Probability Density Plots'
os.makedirs(output_directory, exist_ok=True)

# Iterate through each file in the directory
for file_name in file_list:
    if file_name.endswith('.csv'):
        # Load the viral load data from the Excel file
        data = pd.read_csv(os.path.join(data_directory, file_name))

        # Prepare data for plotting
        viral_loads = np.round(data.values * 10)
        time_steps = np.arange(viral_loads.shape[1])

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
            bin_heights[:, i] = hist / len(viral_loads)  # Compute normalized probability density

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
        # plt.show()

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

        # Save the plots to the output directory & Modify the file name as needed
        plot_file_name = f"{file_name.split('.')[0]}_plot.png"
        plot_path = os.path.join(output_directory, plot_file_name)

        # Save the 3D plot
        fig.savefig(plot_path.replace('.png', '_3D.png'))
        plt.close(fig)

        # Save the 2D plot
        fig2.savefig(plot_path.replace('.png', '_2D.png'))
        plt.close(fig2)


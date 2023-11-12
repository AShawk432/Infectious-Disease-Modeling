import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import warnings

age_groups = ['Overall population', '0-4', '5-14', '15-19', '20-39', '40-59', '60-69', '70-100']

def analyze_csv_directory(directory_path):
    # Create empty lists to store mean and standard deviation for each file
    file_means = []
    file_std_devs = []

    # Iterate through each file in the directory
    for file_name in os.listdir(directory_path):
        if file_name.endswith(".csv"):
            file_path = os.path.join(directory_path, file_name)

            # Read the CSV file into a pandas DataFrame
            df = pd.read_csv(file_path, header=None)

            # Calculate mean and standard deviation for each file
            file_mean = np.mean(df, axis=0).tolist()
            file_std_dev = np.std(df, axis=0).tolist()

            # Append results to lists
            file_means.append(file_mean)
            file_std_devs.append(file_std_dev)

            # Print results for the current file
            print(f"File: {file_name}")
            print("Overall Means:")
            print(file_mean)
            print("\nOverall Standard Deviations:")
            print(file_std_dev)
            print("\n")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        growth_rate = np.divide(file_means[file_index2], file_means[file_index1])

    # Plot the growth rate
    plt.plot(growth_rate, label=f'{age_groups[file_index1]} to {age_groups[file_index2]}', marker='o')

    plt.xlabel('Vector Index')
    plt.ylabel('Growth Rate')
    plt.legend()
    plt.title(f'Overall Mean Growth Rate between age groups ({age_groups[file_index1]} to {age_groups[file_index2]})')
    plt.show()

    # Return the results
    return file_means, file_std_devs

# Read in path to CSV files directory
directory_path =  r'C:\Users\antho\PycharmProjects\pythonProject\Simulation_stat_analysis_data'
file_index1 = 6  # Choose the index of the first file
file_index2 = 2  # Choose the index of the second file

means, std_devs = analyze_csv_directory(directory_path)

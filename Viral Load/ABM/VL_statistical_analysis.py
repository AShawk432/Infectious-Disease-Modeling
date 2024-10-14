import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def plot_variance_and_ci(data, age_group, save_dir=None):
    mean = np.mean(data, axis=0)
    variance = np.var(data, axis=0)
    std_dev = np.std(data, axis=0)
    ci = 1.645 * (std_dev / np.sqrt(len(data)))  # 90% confidence interval using Z-score for normal distribution
    time_steps = np.arange(len(variance))

    plt.plot(time_steps, mean, label='Mean')
    plt.fill_between(time_steps, mean - variance, mean + variance, color='gray', alpha=0.2, label='90% Interval')
    plt.fill_between(time_steps, mean - ci, mean + ci, color='orange', alpha=0.5, label='Variance')
    plt.xlabel('Time Steps')
    plt.ylabel('Viral Load')
    plt.title(f'Viral Load Variance and 90% CI for Age Group: {age_group}')
    plt.legend()
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, f'{age_group}_viral_load_variance_ci.png'))
    else:
        plt.show()
    plt.close()

    return time_steps, ci

def plot_ci_widths(all_ci_data, save_dir=None):
    plt.figure(figsize=(10, 6))
    for age_group, (time_steps, ci_widths) in all_ci_data.items():
        plt.plot(time_steps, ci_widths, label=age_group)

    plt.xlabel('Time Steps')
    plt.ylabel('CI Width')
    plt.title('CI Width for All Age Groups')
    plt.legend()
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, 'all_age_groups_ci_widths.png'))
    else:
        plt.show()
    plt.close()

def plot_ci_ratio(ci_data, age_group1, age_group2, save_dir=None):
    time_steps, ci1 = ci_data[age_group1]
    _, ci2 = ci_data[age_group2]

    # Ensure no division by zero
    ratio = np.divide(ci1, ci2, out=np.zeros_like(ci1), where=ci2!=0)
    slope = np.gradient(ratio, time_steps)  # Calculate the slope of the CI ratio

    plt.figure(figsize=(10, 6))
    plt.plot(time_steps, ratio, label=f'CI Ratio 70-100 and 15-19')
    plt.plot(time_steps, slope, label=f'Slope of CI Ratio 70-100 and 15-19', linestyle='--')
    plt.xlabel('Time Steps')
    plt.ylabel('CI Width Ratio / Slope')
    plt.title(f'CI Width Ratio and Slope Between 70-100 and 15-19')
    plt.legend()
    if save_dir:
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        plt.savefig(os.path.join(save_dir, f'ci_ratio_slope_{age_group1}_vs_{age_group2}.png'))
    else:
        plt.show()
    plt.close()

def main():
    directory = 'C:/Users/antho/PycharmProjects/pythonProject/Primary ABM Model Directory/Viral_Load_Data'  # Update this with the path to your directory of CSV files
    save_directory = 'C:/Users/antho/PycharmProjects/pythonProject/Primary ABM Model Directory/ABM_VL_Plotting'  # Update this with the path where you want to save the plots
    files = os.listdir(directory)
    all_ci_data = {}

    for file in files:
        if file.endswith('.csv'):
            age_group = os.path.splitext(file)[0]
            file_path = os.path.join(directory, file)
            df = pd.read_csv(file_path)

            # Assuming each row represents a person and each column represents a time step
            viral_load_data = df.to_numpy()

            time_steps, ci_widths = plot_variance_and_ci(viral_load_data, age_group, save_directory)
            all_ci_data[age_group] = (time_steps, ci_widths)

    plot_ci_widths(all_ci_data, save_directory)

    # Choose age groups to compare
    age_group1 = 'viral_load_data_by_age_and_time_70-100'  # Replace with actual age group
    age_group2 = 'viral_load_data_by_age_and_time_15-19'  # Replace with actual age group

    if age_group1 in all_ci_data and age_group2 in all_ci_data:
        plot_ci_ratio(all_ci_data, age_group1, age_group2, save_directory)
    else:
        print(f"One or both of the selected age groups ({age_group1}, {age_group2}) are not in the data.")

if __name__ == "__main__":
    main()

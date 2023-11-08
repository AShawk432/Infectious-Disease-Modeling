import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd
import os
import csv

# Define model parameters
num_agents = 500  # Number of agents in the simulation
num_infected = 5  # Number of initially infected agents
# infection_rate = 0.005  # Probability of transmission per contact
latent_period = 5  # Period from getting infected to becoming infectious
infectious_period = 14  # Duration of the infectious period in time steps
time_steps = 60  # Number of time steps in the simulation
immune_period = 7  # Number of days agent is immune from reinfection

# Define age groups and probabilities
age_groups = ['0-4', '5-14', '15-19', '20-39', '40-59', '60-69', '70-100']
age_probs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
# Define death rates by age group
# death_rates = [0.01, 0.02, 0.02, 0.05, 0.1, 0.25, 0.5]
death_rates = [0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07]

# Define the immunosenescence factor for each age group
immunosenescence_factors = [0.95, 0.75, 0.7, 0.5, 0.3, 0.2, 0.1]

# Create a list to store the areas under the viral load curves for each age group
viral_load_areas = []

# Viral load thresholds to determine when agents change compartments
thresh1 = 0.05
thresh2 = 0.9
thresh3 = 0.2

# Define agent class
class Agent:
    def __init__(self, state, viralload, age):
        self.state = state
        self.days_in_compartment = 0
        self.viralload = viralload
        self.immune_days = 0
        self.age = age
        self.is_dead = False
        age_group_index = None
        for index, age_group in enumerate(age_groups):
            age_range = age_group.split('-')
            if int(age_range[0]) <= self.age <= int(age_range[1]):
                age_group_index = index
        self.immunosenescence_factor = immunosenescence_factors[age_group_index]
        self.threshold1 = thresh1 + (random.random() - 0.5) * thresh1
        self.threshold2 = thresh2 + (random.random() - 0.5) * thresh2
        self.threshold3 = thresh3 + (random.random() - 0.5) * thresh3
        self.viral_load_history = []
    def update_state(self, neighbors, deaths_by_ages):
        if self.state == 'S':
            self.days_in_compartment += 1
            # for neighbor in neighbors:
            #     if neighbor.state == 'I' and random.random() < infection_rate:
            #         self.viralload += random.random() / 3
            if self.viralload > self.threshold1:
                self.state = 'E'
                self.days_in_compartment = 0
        elif self.state == 'E':
            self.days_in_compartment += 1
            self.viralload += random.random() / 3
            if self.days_in_compartment < latent_period and self.viralload > self.threshold2:
                self.state = 'I'
                self.days_in_compartment = 0
            elif self.days_in_compartment >= latent_period:
                self.state = 'R'
                self.days_in_compartment = 0

        elif self.state == 'I':
            self.days_in_compartment += 1
            self.viralload -= random.random() * self.immunosenescence_factor
            self.viralload = max(self.viralload, 0)  # Prevent viral load from going below zero
            # Check if agent should die based on age and death rate
            age_group_index = None
            for index, age_group in enumerate(age_groups):
                age_range = age_group.split('-')
                if int(age_range[0]) <= self.age <= int(age_range[1]):
                    age_group_index = index
            if age_group_index is not None:
                if random.random() < death_rates[age_group_index]:
                    self.is_dead = True

            if self.is_dead:
                self.state = 'D'
                self.days_in_compartment = 0
                # Increment deaths in the corresponding age group
                deaths_by_ages[age_group_index] += 1
            if self.viralload <= self.threshold3:
                self.state = 'R'
                self.days_in_compartment = 0

        elif self.state == 'D':
             self.viralload = 0
        elif self.state == 'R':
            self.days_in_compartment += 1
            self.viralload -= self.viralload
            # # ## Adds reinfectivity
            # if self.immune_days >= immune_period:  # Check if the agent's immunity period is over
            #     self.state = 'S'
            #     self.days_in_compartment = 0
            #     self.immune_days = 0    # Reset the immune days counter
            # else:
            #     self.immune_days += 1
            # Append the viral load to the history if it's nonzero
        if self.viralload > 0:
            self.viral_load_history.append(self.viralload)

    def get_state(self):
        return self.state
    def get_age(self):
        return self.age
    def get_age_group(self):
        for age_group in age_groups:
            age_range = age_group.split('-')
            if int(age_range[0]) <= self.age <= int(age_range[1]):
                return age_group
    def die(self):
        self.is_dead = True


# Define simulation function
def simulate():

    # Initialize agents
    agents = []
    people_count = [0] * len(age_groups)
    deaths_by_ages = [0] * len(death_rates)


    # Empty list to append the average viral loads at each time step
    avg_viral_loads = []
    # Initialize a list to store viral load data for each agent at each time step
    viral_load_data_by_agent = []
    for i in range(num_agents):
        if i < num_infected:
            state = 'I'
            viralload = (thresh2 + thresh3) / 2
        else:
            state = 'S'
            viralload = 0

        # Normalize probabilities
        age_probs_normalized = [prob / sum(age_probs) for prob in age_probs]
        # Assign age based on age groups and probabilities
        age_group = np.random.choice(age_groups, p=age_probs_normalized)
        age_range = age_group.split('-')
        age = random.randint(int(age_range[0]), int(age_range[1]))

        agent = Agent(state, viralload, age)
        agents.append(agent)
        viral_load_data_by_agent.append([])
        # Increment the people count for the corresponding age group
        people_count[age_groups.index(age_group)] += 1

    # Run simulation
    state_counts = []
    state_dynamics_by_age = {age_group: [] for age_group in age_groups}  # Dictionary of state dynamics in each age group
    viral_load_data = [[] for _ in range(num_agents)]
    # Create lists to store viral load data for each age group
    viral_load_data_by_age = [[] for _ in range(len(age_groups))]
    # Create a list to store the average viral loads for each age group at each time step
    avg_viral_loads_by_age = [[] for _ in range(len(age_groups))]
    # Create lists to store maximum viral loads for each age group
    max_viral_loads_by_age = [0.0] * len(age_groups)
    # Create a list to store viral load data for each time step and each age group
    viral_load_data_by_age_and_time = [[[] for _ in range(time_steps)] for _ in range(len(age_groups))]
    std_dev_max_viral_loads_by_age = []
    for t in range(time_steps):
        # Update agent states
        for agent in agents:
            neighbors = [neighbor for neighbor in agents if neighbor != agent]
            agent.update_state(neighbors, deaths_by_ages)

            # Get the age group of the current agent
            age_group_index = None
            for index, age_group in enumerate(age_groups):
                age_range = age_group.split('-')
                if int(age_range[0]) <= agent.age <= int(age_range[1]):
                    age_group_index = index

            if age_group_index is not None:
                # Append viral load data to corresponding age group list
                viral_load_data_by_age[age_group_index].append(agent.viralload)
                # print(viral_load_data_by_age)
                # Calculate the maximum viral load for each agent within their age group
        for agent in agents:
            age_group_index = age_groups.index(agent.get_age_group())
            max_viral_loads_by_age[age_group_index] = max(max_viral_loads_by_age[age_group_index], agent.viralload)

        # Create a social interaction matrix based on age group
        social_interaction_matrix = np.array([
            [2.5982, 0.8003, 0.3160, 0.7934, 0.3557, 0.1548, 0.0564],
            [0.6473, 4.1960, 0.6603, 0.5901, 0.4665, 0.1238, 0.0515],
            [0.1737, 1.7500, 11.1061, 0.9782, 0.7263, 0.0815, 0.0273],
            [0.5504, 0.5906, 1.2004, 1.8813, 0.9165, 0.1370, 0.0397],
            [0.3894, 0.7848, 1.3139, 1.1414, 1.3347, 0.2260, 0.0692],
            [0.3610, 0.3918, 0.3738, 0.5248, 0.5140, 0.7072, 0.1469],
            [0.1588, 0.3367, 0.3406, 0.2286, 0.3637, 0.3392, 0.3868]
        ])
        # Normalize the social interaction matrix and compute rolling sums of the rows
        normalized_matrix = social_interaction_matrix / np.sum(social_interaction_matrix, axis=1, keepdims=True)
        row_sums = np.cumsum(normalized_matrix, axis=1)

        # Modify the interaction loop inside the simulation
        for _ in range(500):
            # print("random interaction")
            agent1 = random.choice(agents)  # Choose a random agent

            # Choose the second agent based on age group using the rolling sums
            random_value = random.random()
            age_group1 = agent1.get_age_group()  # Use get_age_group() method
            age_group_index1 = age_groups.index(age_group1)

            # agents_in_age_group = [agent for agent in agents if agent.get_age_group() == age_group]
            probabilities = row_sums[age_group_index1]
            age_group_index2 = np.argmax(
                probabilities > random_value)  # Find the first index where probability exceeds random_value
            age_group2 = age_groups[age_group_index2]  # Get the age group based on the index
            agents_in_age_group2 = [agent for agent in agents if agent.get_age_group() == age_group2]
            if age_group_index2 < len(agents_in_age_group2):
                agent2 = random.choice(agents_in_age_group2)

                # Check if one agent is susceptible and the other is infected
                if (agent1.get_state() == 'S' or agent1.get_state() == 'E') and agent2.get_state() == 'I':
                    susceptible_exposed_agent = agent1
                    infected_agent = agent2
                elif agent1.get_state() == 'I' and (agent2.get_state() == 'S' or agent2.get_state() == 'E'):
                    susceptible_exposed_agent = agent2
                    infected_agent = agent1
                else:
                    continue

                susceptible_exposed_agent.viralload += infected_agent.viralload / 3

        # Record state counts
        s_count = sum([1 for agent in agents if agent.get_state() == 'S'])
        e_count = sum([1 for agent in agents if agent.get_state() == 'E'])
        i_count = sum([1 for agent in agents if agent.get_state() == 'I'])
        r_count = sum([1 for agent in agents if agent.get_state() == 'R'])
        d_count = sum([1 for agent in agents if agent.get_state() == 'D'])
        state_counts.append([s_count, e_count, i_count, r_count, d_count])

        # Calculate state dynamics for each age group
        for age_group in age_groups:
            s_count_age = sum(1 for agent in agents if agent.get_age_group() == age_group and agent.get_state() == 'S')
            e_count_age = sum(1 for agent in agents if agent.get_age_group() == age_group and agent.get_state() == 'E')
            i_count_age = sum(1 for agent in agents if agent.get_age_group() == age_group and agent.get_state() == 'I')
            r_count_age = sum(1 for agent in agents if agent.get_age_group() == age_group and agent.get_state() == 'R')
            d_count_age = sum(1 for agent in agents if agent.get_age_group() == age_group and agent.get_state() == 'D')
            state_dynamics_by_age[age_group].append((s_count_age, e_count_age, i_count_age, r_count_age, d_count_age))

        ## Calculate the average viral load for all agents
        avg_viral_load = sum(agent.viralload for agent in agents if agent.get_state() != 'D') \
            / len([agent for agent in agents if agent.get_state() != 'D'])
        avg_viral_loads.append(avg_viral_load)
        # print(avg_viral_loads)
        # Calculate the standard deviation of the maximum viral loads across all age groups
        std_dev_max_viral_loads_by_age = np.std(max_viral_loads_by_age)

        # Calculate average viral loads for each age group
        for age_group_index, age_group in enumerate(age_groups):
            agents_in_age_group = [agent for agent in agents if
                                   agent.get_age_group() == age_group and agent.get_state() != 'D']
            if agents_in_age_group:
                avg_load_at_time_step = sum(agent.viralload for agent in agents_in_age_group) / len(agents_in_age_group)
                # Update the maximum viral load for the age group
                # avg_max_viral_loads_by_age[age_group_index] = max(avg_max_viral_loads_by_age[age_group_index],
                #                                               avg_load_at_time_step)
            else:
                avg_load_at_time_step = 0  # Handle the case where there are no agents in the age group
            avg_viral_loads_by_age[age_group_index].append(avg_load_at_time_step)


        # Append viral load data for each agent at the current time step
        for i, agent in enumerate(agents):
            viral_load_data[i].append(agent.viralload)
            viral_load_data_by_agent[i].append(agent.viralload)
            age_group_index = age_groups.index(agent.get_age_group())
            viral_load_data_by_age_and_time[age_group_index][t].append(agent.viralload)

    age_df = pd.DataFrame({'Age Group': age_groups, 'People': people_count, 'Deaths': deaths_by_ages})
    print(age_df)

    # Calculate areas under the viral load curves for each age group
    for age_viral_loads in viral_load_data_by_age:
        area_under_curve = np.trapz(age_viral_loads)
        viral_load_areas.append(area_under_curve)

    # Print the areas under the viral load curves for each age group and max avg viral load
    # print("Areas under viral load curves:", viral_load_areas)
    # Print the maximum viral load for each age group
    for age_group_index, age_group in enumerate(age_groups):
        max_viral_load = max_viral_loads_by_age[age_group_index]
        print(f"Maximum Viral Load for {age_group}: {max_viral_load}")
    print(f"Standard Deviation of Maximum Viral Loads: {std_dev_max_viral_loads_by_age}")


    #     # Print ages of all agents
    # for i, agent in enumerate(agents):
    #     print(f"Agent {i + 1} age: {agent.get_age()}")

    # Create a directory to store age group-specific data
    viral_load_dir = "Viral_Load_Data"
    if not os.path.exists(viral_load_dir):
        os.mkdir(viral_load_dir)
    # Write viral load data to a file
    with open(os.path.join(viral_load_dir, 'viral_load.csv'), 'w') as file:
        for agent_loads in viral_load_data:
            file.write(','.join(str(load) for load in agent_loads) + '\n')

    # Modify the code to create separate CSV files for each age group
    for age_group_index, age_group in enumerate(age_groups):
        age_group_file_path = os.path.join(viral_load_dir, f'viral_load_age_{age_group}.csv')
        age_group_data = viral_load_data_by_age_and_time[age_group_index]
        # Transpose the data for this age group
        transposed_data = list(map(list, zip(*age_group_data)))
        with open(age_group_file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            for i, agent_data in enumerate(transposed_data):
                writer.writerow(agent_data)  # Write agent ID and viral load data

    return state_counts, agents, avg_viral_loads, state_dynamics_by_age, avg_viral_loads_by_age, viral_load_data_by_age, \



# # Run simulation
# state_counts, agents, avg_viral_loads, viral_load_data_by_agent = simulate()
# state_counts = np.array(state_counts)
#
# s_counts = state_counts[:, 0]
# e_counts = state_counts[:, 1]
# i_counts = state_counts[:, 2]
# r_counts = state_counts[:, 3]
# d_counts = state_counts[:, 4]

# Run simulation 100 times and accumulate results
num_simulations = 2
avg_state_counts = np.zeros((time_steps, 5))  # Initialize an array to accumulate state counts
overall_avg_loads = []
avg_state_dynamics_by_age = {age_group: [] for age_group in age_groups}
avg_viral_load_by_age = [[] for _ in range(len(age_groups))]
overall_avg_loads_by_age = []
simulation_data_by_age_group = {age_group: [] for age_group in age_groups}
viral_load_histories_by_age_group = {age_group: [] for age_group in age_groups}
viral_load_histories = []

for simulation in range(num_simulations):

    # print(avg_viral_loads)
    state_counts, agents, avg_viral_loads, state_dynamics_by_age, avg_viral_loads_by_age, viral_load_data_by_age, \
         = simulate()

    avg_state_counts += np.array(state_counts)
    # Store the average viral loads and profiles at each time step for this simulation
    overall_avg_loads.append(avg_viral_loads)
    for agent in agents:
        age_group = agent.get_age_group()
        viral_load_history = agent.viral_load_history
        viral_load_histories_by_age_group[age_group].append(viral_load_history)
        max_history_length = max(len(history) for history in viral_load_histories_by_age_group.values())
        for age_group in viral_load_histories_by_age_group:
            viral_load_histories_by_age_group[age_group] += [0] * (
                        max_history_length - len(viral_load_histories_by_age_group[age_group]))
    # padded_viral_load_histories = []
    # for age_group, viral_load_history in viral_load_histories_by_age_group.items():
    #     padding_length = max_history_length - len(viral_load_history)
    #     padded_viral_load_history = viral_load_history + [0] * padding_length
    #     padded_viral_load_histories.append(padded_viral_load_history)


    for age_group in age_groups:
        age_group_index = age_groups.index(age_group)
        overall_avg_loads_by_age.append(avg_viral_loads_by_age)
        simulation_data_by_age_group[age_group].append(avg_viral_loads_by_age[age_group_index])

    # Accumulate state dynamics by age
    for age_group in age_groups:
        avg_state_dynamics_by_age[age_group].append(np.array(state_dynamics_by_age[age_group]))

# Create a directory to store age group-specific data
viral_load_data_dir = "Simulation_stat_analysis_data"
if not os.path.exists(viral_load_data_dir):
    os.mkdir(viral_load_data_dir)
# Save the overall average viral load data to a CSV file in the same directory as age group data
overall_avg_file_path = os.path.join(viral_load_data_dir, "overall_avg_viral_load.csv")
with open(overall_avg_file_path, 'w', newline='') as overall_file:
    writer = csv.writer(overall_file)
    writer.writerows(overall_avg_loads)
# Write the data for each age group to separate CSV files
for age_group_index, age_group in enumerate(age_groups):
    age_group_file_path = os.path.join(viral_load_data_dir, f'overall_viral_load_age_{age_group}.csv')
    age_group_data = np.array(simulation_data_by_age_group[age_group])
    with open(age_group_file_path, 'w', newline='') as age_file:
        writer = csv.writer(age_file)
        # header_row = [str(i) for i in range(age_group_data.shape[1])]
        # writer.writerow(header_row)
        writer.writerows(age_group_data)

# Calculate the overall average viral load at each time step across all simulations
overall_avg_viral_loads = np.mean(np.array(overall_avg_loads), axis=0)
for age_group in age_groups:
    overall_avg_viral_loads_by_age = np.mean(np.array(overall_avg_loads_by_age), axis=0)
    avg_viral_load_profiles_by_age_group = np.mean(np.array(list(viral_load_histories_by_age_group.values())), axis=0)

# Calculate the average state dynamics by age
for age_group in age_groups:
    avg_state_dynamics_by_age[age_group] = np.mean(np.array(avg_state_dynamics_by_age[age_group]), axis=0)

# Extract individual state counts for plotting
s_counts = avg_state_counts[:, 0]
e_counts = avg_state_counts[:, 1]
i_counts = avg_state_counts[:, 2]
r_counts = avg_state_counts[:, 3]
d_counts = avg_state_counts[:, 4]


def plotting_function():

    # Plot SEIR dynamics for each state of agents over time
    plt.figure(figsize=(10, 8))
    plt.plot(s_counts, label='Susceptible')
    plt.plot(e_counts, label='Exposed')
    plt.plot(i_counts, label='Infected')
    plt.plot(r_counts, label='Recovered')
    plt.plot(d_counts, label='Deaths')
    plt.xlabel('Time steps')
    plt.ylabel('Number of agents')
    plt.title('Agent-based SEIRD model simulation')
    plt.legend()
    plt.grid(True)
    plt.show()

    # Collect time total steps in a vector
    step_count = []
    for steps in range(time_steps):
        step_count.append(steps)

    # Plot the average viral loads over time
    plt.figure(figsize=(10, 8))
    plt.plot(step_count, avg_viral_loads, label='Total Viral Load', color='purple')
    plt.title('Average Viral Load Over Time')
    plt.xlabel('Time Steps')
    plt.ylabel('Viral Load')
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.legend()
    plt.grid(True)
    plt.show()

    # Plot the average viral loads over time
    plt.figure(figsize=(10, 8))
    plt.plot(step_count, overall_avg_viral_loads, label='Average Viral Load', color='purple')
    plt.title('Average Viral Load Over Time (Averaged Across Simulations)')
    plt.xlabel('Time Steps')
    plt.ylabel('Average Viral Load')
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.legend()
    plt.grid(True)
    plt.show()

    # Create a directory to store age group state dynamics plots
    plotting_dir = "ABM_VL_Plotting"
    if not os.path.exists(plotting_dir):
        os.mkdir(plotting_dir)

    # Plot state dynamics for each age group and save to the folder
    for age_group in age_groups:
        dynamics_data = avg_state_dynamics_by_age[age_group]
        s_counts_age = [data[0] for data in dynamics_data]
        e_counts_age = [data[1] for data in dynamics_data]
        i_counts_age = [data[2] for data in dynamics_data]
        r_counts_age = [data[3] for data in dynamics_data]
        d_counts_age = [data[4] for data in dynamics_data]

        plt.figure(figsize=(10, 8))
        plt.plot(s_counts_age, label='Susceptible')
        plt.plot(e_counts_age, label='Exposed')
        plt.plot(i_counts_age, label='Infected')
        plt.plot(r_counts_age, label='Recovered')
        plt.plot(d_counts_age, label='Deaths')
        plt.xlabel('Time steps')
        plt.ylabel('Number of agents')
        plt.title(f'State Dynamics for Age Group {age_group}')
        plt.legend()
        plt.grid(True)

        # Save the state dynamics plots for each age group into the folder
        plt.savefig(os.path.join(plotting_dir, f'age_group_{age_group}_step_{time_steps}.png'))
        plt.close()

    for age_group_index, age_group in enumerate(age_groups):
        plt.figure(figsize=(10, 8))
        plt.plot(overall_avg_viral_loads_by_age[age_group_index], label=f'Age Group {age_group}', color='red')
        plt.xlabel('Time steps')
        plt.ylabel('Average Viral Load')
        plt.title(f'Average Viral Load for Age Group {age_group} Over Time')
        plt.legend()
        plt.grid(True)
        avg_viral_loads_filename = f'average_viral_loads_age_group_{age_group}.png'
        avg_viral_loads_filepath = os.path.join(plotting_dir, avg_viral_loads_filename)
        plt.savefig(avg_viral_loads_filepath)
        plt.close()

        # Plot the viral load curves for each age group on the same plot with different colors
    plt.figure(figsize=(10, 8))
    for age_group_index, age_group in enumerate(age_groups):
        plt.plot(step_count, overall_avg_viral_loads_by_age[age_group_index], label=f'Age Group {age_group}', alpha=0.7)
    plt.xlabel('Time steps')
    plt.ylabel('Average Viral Load')
    plt.title('Average Viral Load Over Time by Age Group')
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.legend()
    plt.grid(True)
    plt.show()

    # Calculate and plot the average viral load profile for each age group
    for age_group in age_groups:
        avg_viral_load_profile = viral_load_histories_by_age_group[age_group]
        plt.figure(figsize=(10, 8))
        plt.plot(avg_viral_load_profile, label=f'Age Group {age_group}', color='green')
        plt.xlabel('Time steps')
        plt.ylabel('Average Viral Load Profile')
        plt.title(f'Average Viral Load Profile for Age Group {age_group}')
        plt.legend()
        plt.grid(True)
        plt.savefig(os.path.join(plotting_dir, f'viral_load_profile_age_group_{age_group}.png'))
        plt.close()

        plt.plot(avg_viral_load_profile, label=f'Age Group {age_group}')
    plt.xlabel('Time steps')
    plt.ylabel('Average Viral Load Profile')
    plt.title('Average Viral Load Profiles for All Age Groups')
    plt.legend()
    plt.grid(True)
    plt.savefig(os.path.join(plotting_dir, f'Average Viral Load Profiles for All Age Groups.png'))


plotting_function()

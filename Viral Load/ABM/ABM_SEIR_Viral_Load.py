import numpy as np
import random
import matplotlib.pyplot as plt
import pandas as pd

# Define model parameters
num_agents = 500  # Number of agents in the simulation
num_infected = 5  # Number of initially infected agents
# infection_rate = 0.005  # Probability of transmission per contact
latent_period = 5  # Period from getting infected to becoming infectious
infectious_period = 14  # Duration of the infectious period in time steps
time_steps = 100  # Number of time steps in the simulation
immune_period = 7  # Number of days agent is immune from reinfection

# Define age groups and probabilities
age_groups = ['0-4', '5-14', '15-19', '20-39', '40-59', '60-69', '70-100']
age_probs = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]
# Define death rates by age group
# death_rates = [0.01, 0.02, 0.02, 0.05, 0.1, 0.25, 0.5]
death_rates = [0.07, 0.07, 0.07, 0.07, 0.07, 0.07, 0.07]

# Define the immunosenescence factor for each age group
immunosenescence_factors = [0.95, 0.75, 0.7, 0.5, 0.3, 0.2, 0.1]

# Empty list to append the average viral loads at each time step
viral_loads = []
# Create a list to store the areas under the viral load curves for each age group
viral_load_areas = []

# Viral load thresholds to determine when agents change compartments
thresh1 = 0.05
thresh2 = 0.9
thresh3 = 0.2

matrix = [0,0,0,0,1,3]
index = np.argmax(matrix)
print(index)
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
    def update_state(self, neighbors, deaths_by_ages):
        if self.state == 'S':
            self.days_in_compartment += 1
            # for neighbor in neighbors:
            #     if neighbor.state == 'I' and random.random() < infection_rate:
            #         self.viralload += random.random() / 3
            if self.viralload > thresh1:
                self.state = 'E'
                self.days_in_compartment = 0
        elif self.state == 'E':
            self.days_in_compartment += 1
            self.viralload += random.random() / 3
            if self.days_in_compartment < latent_period and self.viralload > thresh2:
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
            if self.viralload <= thresh3:
                self.state = 'R'
                self.days_in_compartment = 0

        elif self.state == 'D':
             self.viralload = 0
        elif self.state == 'R':
            self.days_in_compartment += 1
            self.viralload -= self.viralload
            # # ## Adds reinfectivity
            # if self.viralload <= 0:
            #     if self.immune_days >= immune_period:  # Check if the agent's immunity period is over
            #         self.state = 'S'
            #         self.days_in_compartment = 0
            #         self.immune_days = 0    # Reset the immune days counter
            #     else:
            #         self.immune_days += 1

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
    # Create lists to store viral load data for each age group
    viral_load_data_by_age = [[] for _ in range(len(age_groups))]
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

        # Increment the people count for the corresponding age group
        people_count[age_groups.index(age_group)] += 1

    # Run simulation
    state_counts = []
    viral_load_data = [[] for _ in range(num_agents)]
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

        # social_interaction_matrix = np.array([
        #     [0.1, 0.05, 0.03, 0.01, 0.01, 0.005, 0.001],
        #     [0.05, 0.1, 0.05, 0.03, 0.02, 0.01, 0.005],
        #     [0.03, 0.05, 0.1, 0.05, 0.03, 0.02, 0.01],
        #     [0.01, 0.03, 0.05, 0.1, 0.05, 0.03, 0.02],
        #     [0.01, 0.02, 0.03, 0.05, 0.1, 0.05, 0.03],
        #     [0.005, 0.01, 0.02, 0.03, 0.05, 0.1, 0.05],
        #     [0.001, 0.005, 0.01, 0.02, 0.03, 0.05, 0.1]
        # ]

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
                if agent1.get_state() == 'S' and agent2.get_state() == 'I':
                    susceptible_agent = agent1
                    infected_agent = agent2
                elif agent1.get_state() == 'I' and agent2.get_state() == 'S':
                    susceptible_agent = agent2
                    infected_agent = agent1
                else:
                    continue

                susceptible_agent.viralload += infected_agent.viralload / 3

        # Record state counts
        s_count = sum([1 for agent in agents if agent.get_state() == 'S'])
        e_count = sum([1 for agent in agents if agent.get_state() == 'E'])
        i_count = sum([1 for agent in agents if agent.get_state() == 'I'])
        r_count = sum([1 for agent in agents if agent.get_state() == 'R'])
        d_count = sum([1 for agent in agents if agent.get_state() == 'D'])
        state_counts.append([s_count, e_count, i_count, r_count, d_count])


        ## Calculate the average viral load for all agents
        avg_viral_loads = sum(agent.viralload for agent in agents if agent.get_state() != 'D') \
            / len([agent for agent in agents if agent.get_state() != 'D'])
        viral_loads.append(avg_viral_loads)
        # print(avg_viral_loads)

        # Append viral load data for each agent at the current time step
        for i, agent in enumerate(agents):
            viral_load_data[i].append(agent.viralload)

    age_df = pd.DataFrame({'Age Group': age_groups, 'People': people_count, 'Deaths': deaths_by_ages})
    print(age_df)

    # Calculate areas under the viral load curves for each age group
    for age_viral_loads in viral_load_data_by_age:
        area_under_curve = np.trapz(age_viral_loads)
        viral_load_areas.append(area_under_curve)

    # Print the areas under the viral load curves for each age group
    print("Areas under viral load curves:", viral_load_areas)

    #     # Print ages of all agents
    # for i, agent in enumerate(agents):
    #     print(f"Agent {i + 1} age: {agent.get_age()}")

    # # Write viral load data to a file
    # with open('viral_load.csv', 'w') as file:
    #     for agent_loads in viral_load_data:
    #         file.write(','.join(str(load) for load in agent_loads) + '\n')

    return state_counts, agents, viral_loads

# Run simulation
state_counts, agents, viral_loads = simulate()
state_counts = np.array(state_counts)

s_counts = state_counts[:, 0]
e_counts = state_counts[:, 1]
i_counts = state_counts[:, 2]
r_counts = state_counts[:, 3]
d_counts = state_counts[:, 4]


# # Run simulation 100 times and accumulate results
# num_simulations = 2
# avg_state_counts = np.zeros((time_steps, 5))  # Initialize an array to accumulate state counts
# viral_loads_matrix = np.zeros((num_simulations, time_steps))  # Initialize a matrix for viral load data
#
# for simulation in range(num_simulations):
#     state_counts, agents, viral_loads = simulate()
#     avg_state_counts += np.array(state_counts)
#     viral_loads_matrix[simulation, :] = viral_loads  # Store the viral load data for this simulation
#
# avg_state_counts /= num_simulations  # Calculate the average state counts over all simulations
#
# # Calculate the average viral load over simulations for each time step
# avg_viral_loads = np.mean(viral_loads_matrix, axis=0)
#
# # Extract individual state counts for plotting
# s_counts = avg_state_counts[:, 0]
# e_counts = avg_state_counts[:, 1]
# i_counts = avg_state_counts[:, 2]
# r_counts = avg_state_counts[:, 3]
# d_counts = avg_state_counts[:, 4]


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

    ## Collect time total steps in a vector
    step_count = []
    for steps in range(time_steps):
        step_count.append(steps)

    # Plot the average viral loads over time
    plt.figure(figsize=(10, 8))
    plt.plot(step_count, viral_loads, label='Total Viral Load', color='purple')
    plt.title('Average Viral Load Over Time')
    plt.xlabel('Time Steps')
    plt.ylabel('Viral Load')
    plt.xticks(rotation=45)
    plt.yticks(rotation=45)
    plt.legend()
    plt.grid(True)
    plt.show()

# Call the new function to run the simulation and plot the results
plotting_function()

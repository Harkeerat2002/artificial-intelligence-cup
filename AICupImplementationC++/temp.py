import numpy as np
import random

# Read in the .tsp file and extract the necessary information
with open('eil76.tsp', 'r') as f:
    lines = f.readlines()
    dimension = int(lines[3].split()[1])
    coordinates = []
    for i in range(6, len(lines)):
        x, y = map(float, lines[i].split())
        coordinates.append((x, y))

# Calculate the Euclidean distance between two points
def euclidean_distance(a, b):
    return np.sqrt((a[0] - b[0])**2 + (a[1] - b[1])**2)

# Calculate the distances between all pairs of cities
distances = np.zeros((dimension, dimension))
for i in range(dimension):
    for j in range(dimension):
        distances[i][j] = euclidean_distance(coordinates[i], coordinates[j])

# Parameters for the ACO algorithm
num_ants = 50
num_iterations = 100
alpha = 1
beta = 5
rho = 0.1

# Initialize the pheromone matrix
pheromones = np.ones((dimension, dimension)) / dimension

# Initialize the heuristic information matrix
heuristics = 1 / distances

# Initialize the current solutions for each ant
current_solutions = np.zeros((num_ants, dimension))
current_solutions[:, 0] = np.arange(num_ants)

# Initialize the current costs for each ant
current_costs = np.zeros(num_ants)

# Initialize the best solution and cost
best_solution = np.zeros(dimension)
best_cost = float('inf')

# ACO loop
for it in range(num_iterations):
    # Update the current solutions for each ant
    for a in range(num_ants):
        for i in range(1, dimension):
            # Calculate the probability of choosing each city
            prob = (pheromones[int(current_solutions[a][i-1])][:] ** alpha) * (heuristics[int(current_solutions[a][i-1])][:] ** beta)
            prob /= sum(prob)
            # Choose the next city based on the probabilities
            current_solutions[a][i] = np.random.choice(dimension, p=prob)
        # Calculate the cost of the current solution for this ant
        current_costs[a] = sum([distances[int(current_solutions[a][i])][int(current_solutions[a][i+1])] for i in range(dimension-1)])
        current_costs[a] += distances[int(current_solutions[a][dimension-1])][int(current_solutions[a][0])]
        # Update the best solution and cost if necessary
        if current_costs[a] < best_cost:
            best_cost = current_costs[a]
            best_solution = current_solutions[a].copy()
            # Update the pheromones on the edges based on the current solutions
        delta_pheromones = np.zeros((dimension, dimension))

        for a in range(num_ants):
            for i in range(dimension-1):
                delta_pheromones[int(current_solutions[a][i])][int(current_solutions[a][i+1])] += 1 / current_costs[a]
            delta_pheromones[int(current_solutions[a][dimension-1])][int(current_solutions[a][0])] += 1 / current_costs[a]
        pheromones = (1 - rho) * pheromones + delta_pheromones

# Print the best solution and cost for this iteration
print(f'Iteration {it+1}: Best solution = {best_solution}, Best cost = {best_cost}')

# Print the final best solution and cost after all iterations
print(f'Final solution: {best_solution}, Final cost: {best_cost}')


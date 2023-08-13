import numpy as np
from tqdm import tqdm

# Read the .tsp file and store the coordinates of the cities
cities = []
with open('eil76.tsp') as f:
    for line in f:
        if 'NODE_COORD_SECTION' in line:
            break
    for line in f:
        if 'EOF' in line:
            break
        parts = line.split()
        cities.append((float(parts[1]), float(parts[2])))

# Calculate the distance between each pair of cities
num_cities = len(cities)
distances = np.zeros((num_cities, num_cities))
for i in range(num_cities):
    for j in range(i+1, num_cities):
        x1, y1 = cities[i]
        x2, y2 = cities[j]
        distance = np.sqrt((x1-x2)**2 + (y1-y2)**2)
        distances[i][j] = distance
        distances[j][i] = distance

# Initialize the pheromone matrix with a small positive value
pheromones = np.ones((num_cities, num_cities)) * 0.1

# Set the hyperparameters for the ant colony optimization
alpha = 1.0  # pheromone influence
beta = 2.0  # distance influence
evaporation_rate = 0.5  # pheromone evaporation rate
ant_count = 20  # number of ants
max_iterations = 100  # number of iterations

# Create a function to calculate the probability of an ant choosing a particular city as the next move


def calculate_probabilities(ant, visited, pheromones, distances):
    # Extract the current city and the list of unvisited cities for the given ant
    current_city = visited[ant][-1]
    unvisited_cities = [i for i in range(num_cities) if i not in visited[ant]]

    # Initialize an array to store the probabilities for each city
    probabilities = np.zeros(num_cities)

    # Calculate the denominator for the probability formula
    denominator = 0
    for city in unvisited_cities:
        denominator += pheromones[current_city][city]**alpha * \
            (1/distances[current_city][city])**beta

    # Calculate the probability for each unvisited city
    for city in unvisited_cities:
        probabilities[city] = (pheromones[current_city][city]**alpha *
                               (1/distances[current_city][city])**beta) / denominator

    return probabilities

# Create a function to move the ants to the next city


def move_ants(probabilities, visited, pheromones, distances):
    for ant in range(ant_count):
        # Extract the current city and the list of unvisited cities for the given ant
        current_city = visited[ant][-1]


        # Initialize the pheromone update matrix with zeros
pheromone_updates = np.zeros((num_cities, num_cities))

# Run the ant colony optimization algorithm for a specified number of iterations
for iteration in tqdm(range(max_iterations)):
    # Initialize a list to store the cities visited by each ant
    visited = [[] for _ in range(ant_count)]

    # Initialize the ants at the first city
    for ant in range(ant_count):
        visited[ant].append(0)

    # Move the ants to the next city until all cities have been visited
    while len(visited[0]) < num_cities:
        # Calculate the probabilities of each ant moving to a particular city
        probabilities = [calculate_probabilities(
            ant, visited, pheromones, distances) for ant in range(ant_count)]

        # Move the ants to the next city based on the probabilities
        visited = move_ants(probabilities, visited, pheromones, distances)

    # Update the pheromones based on the paths taken by the ants
    for ant in range(ant_count):
        for i in range(num_cities-1):
            pheromone_updates[visited[ant][i]][visited[ant][i+1]
                                               ] += 1.0 / distances[visited[ant][i]][visited[ant][i+1]]
        pheromone_updates[visited[ant][-1]][visited[ant][0]
                                            ] += 1.0 / distances[visited[ant][-1]][visited[ant][0]]

    # Evaporate the pheromones
    pheromones = (1 - evaporation_rate) * pheromones + pheromone_updates
    pheromone_updates = np.zeros((num_cities, num_cities))

# Find the shortest path among all the paths taken by the ants
best_path = min(visited, key=lambda x: sum(
    [distances[x[i]][x[i+1]] for i in range(num_cities-1)]))

# Print the shortest path and the total distance
print(best_path)
print(sum([distances[best_path[i]][best_path[i+1]]
      for i in range(num_cities-1)]))

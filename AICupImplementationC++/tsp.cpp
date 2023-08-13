#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <vector>

// Struct to represent a city in the TSP
struct City
{
    int x, y;
};

// Function to calculate the distance between two cities
double distance(City a, City b)
{
    int dx = a.x - b.x;
    int dy = a.y - b.y;
    return sqrt(dx * dx + dy * dy);
}

// Class to represent an ant in the ACO algorithm
class Ant
{
public:
    // Constructor
    Ant(int num_cities)
    {
        visited.resize(num_cities);
        std::fill(visited.begin(), visited.end(), false);
        path_length = 0;
    }

    // Function to move the ant to a new city
    void moveTo(City city, double **pheromones, double **distances, double alpha, double beta)
    {
        // Mark the city as visited
        visited[city.x] = true;

        // If this is not the first city, update the path length
        if (current_city.x != -1)
        {
            path_length += distance(current_city, city);
        }

        // Update the current city
        current_city = city;

        // Calculate the probabilities of moving to each city
        std::vector<double> probabilities;
        for (int i = 0; i < num_cities; i++)
        {
            if (visited[i])
            {
                probabilities.push_back(0);
            }
            else
            {
                double p = pow(pheromones[current_city.x][i], alpha) * pow(1.0 / distances[current_city.x][i], beta);
                probabilities.push_back(p);
            }
        }

        // Normalize the probabilities
        double sum = 0;
        for (double p : probabilities)
            sum += p;
        for (int i = 0; i < num_cities; i++)
            probabilities[i] /= sum;

        // Choose the next city based on the probabilities
        double r = (double)rand() / RAND_MAX;
        double total = 0;
        for (int i = 0; i < num_cities; i++)
        {
            total += probabilities[i];
            if (total >= r)
            {
                moveTo(cities[i], pheromones, distances, alpha, beta);
                break;
            }
        }
    }

    // Function to reset the ant's state
    void reset()
    {
        std::fill(visited.begin(), visited.end(), false);
        path_length = 0;
        current_city = {-1, -1};
    }

    // Variables to store the state of the ant
    City current_city;
    std::vector<bool> visited;
    double path_length;
    std::vector<City> cities;
    int num_cities;
};

int main()
{
    // Seed the random
    // Seed the random number generator
    srand(time(0));

    // Read in the cities from the .tsp file
    std::vector<City> cities;
    int num_cities;
    std::cin >> num_cities;
    for (int i = 0; i < num_cities; i++)
    {
        int x, y;
        std::cin >> x >> y;
        cities.push_back({x, y});
    }

    // Initialize the pheromones and distances matrices
    double **pheromones = new double *[num_cities];
    double **distances = new double *[num_cities];
    for (int i = 0; i < num_cities; i++)
    {
        double **pheromones = new double *[num_cities];
        for (int i = 0; i < num_cities; i++)
        {
            pheromones[i] = new double[num_cities];
        }

        distances[i] = new double[num_cities];
        for (int j = 0; j < num_cities; j++)
        {
            pheromones[i][j] = 1;
            distances[i][j] = distance(cities[i], cities[j]);
        }
    }

    // Set the parameters for the ACO algorithm
    double alpha = 1;
    double beta = 5;
    int num_ants = 20;
    int num_iterations = 100;
    double rho = 0.1;

    // Create the ants
    std::vector<Ant> ants;
    for (int i = 0; i < num_ants; i++)
    {
        ants.push_back(Ant(num_cities));
        ants[i].cities = cities;
        ants[i].num_cities = num_cities;
    }

    // Run the ACO algorithm
    for (int it = 0; it < num_iterations; it++)
    {
        // Reset the ants
        for (Ant &ant : ants)
            ant.reset();
        // Have each ant construct a solution
        for (Ant &ant : ants)
        {
            ant.moveTo(cities[0], pheromones, distances, alpha, beta);
        }

        // Update the pheromones
        for (int i = 0; i < num_cities; i++)
        {
            for (int j = 0; j < num_cities; j++)
            {
                double delta_tau = 0;
                for (const Ant &ant : ants)
                {
                    if (ant.visited[j])
                    {
                        delta_tau += 1.0 / ant.path_length;
                    }
                }
                pheromones[i][j] *= (1 - rho);
                pheromones[i][j] += delta_tau;
            }
        }
    }

    // Find the best solution
    Ant best_ant = ants[0];
    for (const Ant &ant : ants)
    {
        if (ant.path_length < best_ant.path_length)
        {
            best_ant = ant;
        }
    }

    // Print the best solution
    std::cout << "Best solution found:" << std::endl;
    for (int i = 0; i < num_cities; i++)
    {
        std::cout << best_ant.cities[i].x << " " << best_ant.cities[i].y << std::endl;
    }
    std::cout << "Total distance: " << best_ant.path_length << std::endl;

    // Free the memory
    for (int i = 0; i < num_cities; i++)
    {
        for (int i = 0; i < num_cities; i++)
        {
            delete[] pheromones[i];
        }
        delete[] pheromones;

        for (int i = 0; i < num_cities; i++)
        {
            delete[] distances[i];
        }
        delete[] distances[i];
    }
    delete[] pheromones;
    delete[] distances;

    return 0;
}
#include <iostream>
#include <algorithm>
#include <random>
#include <vector>
#include <cmath>
#include <fstream>
#include <sstream>
#include <map>

const int MAX_N = 130;           // maximum number of cities
const int MAX_ITERATIONS = 1000; // maximum number of iterations
const int INITIAL_PHEROMONE = 1; // initial pheromone value
const int ALPHA = 1;             // importance of pheromone
const int BETA = 1;              // importance of distance
const double RHO = 0.5;          // pheromone decay rate
const double Q = 100;            // pheromone deposit rate

int n;                                           // number of cities
double dist[MAX_N][MAX_N];                       // distance between cities
double pheromone[MAX_N][MAX_N];                  // pheromone on the path between cities
std::map<int, std::pair<double, double>> coords; // coordinates of each city

std::mt19937 rng((unsigned)time(nullptr));         // random number generator
std::uniform_real_distribution<double> unif(0, 1); // random number distribution

// Calculate the distance between two cities
double distance(int i, int j)
{
    auto coord1 = coords[i];
    auto coord2 = coords[j];
    return sqrt((coord1.first - coord2.first) * (coord1.first - coord2.first) + (coord1.second - coord2.second) * (coord1.second - coord2.second));
}

// Calculate the probability of choosing a particular edge
double probability(int i, int j)
{
    double numerator = pow(pheromone[i][j], ALPHA) * pow(1.0 / dist[i][j], BETA);
    double denominator = 0;
    for (int k = 0; k < n; k++)
    {
        if (k != i)
        {
            denominator += pow(pheromone[i][k], ALPHA) * pow(1.0 / dist[i][k], BETA);
        }
    }
    return numerator / denominator;
}

// Find the shortest path using Ant Colony Optimization
std::vector<int> ant_colony_optimization()
{
    // Initialize pheromone on all edges to INITIAL_PHEROMONE
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            pheromone[i][j] = INITIAL_PHEROMONE;
        }
    }

    // Initialize best path to the empty path
    std::vector<int> best_path;
    double best_length = 0;

    // Repeat until maximum number of iterations is reached
    for (int t = 0; t < MAX_ITERATIONS; t++)
    {
        // Print the Iteration
        std::cout << "Iteration: " << t << std::endl;
        // Initialize array of paths
        std::vector<int> path[MAX_N];

        // Initialize array of path lengths
        double path_length[MAX_N];
        // Create paths for each ant
        for (int i = 0; i < n; i++)
        {
            // Choose starting city
            int current_city = i;
            path[i].push_back(current_city);

            // Initialize visited array
            bool visited[MAX_N];
            for (int j = 0; j < n; j++)
            {
                visited[j] = false;
            }
            visited[current_city] = true;

            // Repeat until all cities are visited
            while (path[i].size() < n)
            {
                // Choose next city
                double r = unif(rng);
                double sum = 0;
                for (int j = 0; j < n; j++)
                {
                    if (!visited[j])
                    {
                        sum += probability(current_city, j);
                        if (sum >= r)
                        {
                            current_city = j;
                            break;
                        }
                    }
                }

                // Add next city to path
                path[i].push_back(current_city);
                visited[current_city] = true;
            }

            // Calculate path length
            path_length[i] = 0;
            for (int j = 0; j < n - 1; j++)
            {
                path_length[i] += dist[path[i][j]][path[i][j + 1]];
            }
            path_length[i] += dist[path[i][n - 1]][path[i][0]];

            // Update best path and best length if necessary
            if (best_path.empty() || path_length[i] < best_length)
            {
                best_path = path[i];
                best_length = path_length[i];
            }
        }

        // Deposit pheromone on the best path
        for (int i = 0; i < n - 1; i++)
        {
            pheromone[best_path[i]][best_path[i + 1]] += Q / best_length;
            pheromone[best_path[i + 1]][best_path[i]] += Q / best_length;
        }
        pheromone[best_path[n - 1]][best_path[0]] += Q / best_length;
        pheromone[best_path[0]][best_path[n - 1]] += Q / best_length;

        // Evaporate pheromone on all edges
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                pheromone[i][j] *= 1 - RHO;
            }
        }
    }

    return best_path;
}

int main()
{
    std::ifstream infile("eil76.tsp");
    std::string line;
    int count = 0;
    while (std::getline(infile, line))
    {
        std::cout << line << std::endl;
        std::istringstream iss(line);
        std::string word;
        iss >> word;
        if (word == "DIMENSION:")
        {
            iss >> n;
        }
        else if (word == "NODE_COORD")
        {
            // Read city coordinates
            int city;
            double x, y;
            iss >> city >> x >> y;
            coords[city - 1] = std::make_pair(x, y);

            // Print city coordinates
            std::cout << city << " " << x << " " << y << std::endl;
        }
    }

    std::cout << "Number of cities: " << n << std::endl;

    // Calculate distance between cities
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            dist[i][j] = distance(i, j);
        }
    }

    // print the distance between cities
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << dist[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Find the shortest path using Ant Colony Optimization
    std::vector<int> best_path = ant_colony_optimization();

    // Print the shortest path
    for (int i = 0; i < n; i++)
    {
        std::cout << best_path[i] << " ";
    }
    std::cout << std::endl;

    return 0;
}


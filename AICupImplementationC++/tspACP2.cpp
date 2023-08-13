#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>
#include <algorithm>

using namespace std;

// define the city structure
struct City
{
    int id;
    double x;
    double y;
};

// parse the input data
void parse_input(vector<City>& cities, int& num_cities, double& best_known)
{
    ifstream input_file("ch130.tsp");
    string line;
    bool node_coord_section = false;
    while (getline(input_file, line))
    {
        if (line == "NODE_COORD_SECTION")
        {
            node_coord_section = true;
            continue;
        }
        if (line == "EOF")
        {
            break;
        }

        if (node_coord_section)
        {
            int id;
            double x, y;
            sscanf(line.c_str(), "%d %lf %lf", &id, &x, &y);
            cities.push_back({id, x, y});
        }
        else
        {
            if (line.find("DIMENSION") != string::npos)
            {
                sscanf(line.c_str(), "DIMENSION : %d", &num_cities);
            }
            if (line.find("BEST_KNOWN") != string::npos)
            {
                sscanf(line.c_str(), "BEST_KNOWN : %lf", &best_known);
            }
        }
    }
}

// calculate the Euclidean distance between two cities
double euc_2d_distance(City city1, City city2)
{
    double dx = city1.x - city2.x;
    double dy = city1.y - city2.y;
    return sqrt(dx * dx + dy * dy);
}

// initialize the p
// initialize the pheromone trail matrix
vector<vector<double>> initialize_pheromones(int num_cities)
{
    // set the initial pheromone value to a high value
    double initial_pheromone = 100.0;

    // create a 2D vector with the desired size
    vector<vector<double>> pheromones(num_cities, vector<double>(num_cities, initial_pheromone));

    // set the pheromone for the same city to 0
    for (int i = 0; i < num_cities; i++)
    {
        pheromones[i][i] = 0;
    }

    return pheromones;
}

// find the shortest path using ACO
pair<vector<int>, double> find_shortest_path(vector<City> cities, int num_cities, 
                                              vector<vector<double>> pheromones, int num_ants, 
                                              int num_iterations, double evap_rate)
{
    // create a random number generator
    mt19937 rng(random_device{}());

    // initialize the best solution variables
    vector<int> best_solution;
    double best_distance = numeric_limits<double>::max();

    // run the ACO algorithm for the specified number of iterations
    for (int iteration = 0; iteration < num_iterations; iteration++)
    {
        cout << "Iteration: " << iteration << endl;

        // create a 2D vector to store the paths taken by each ant
        vector<vector<int>> paths(num_ants, vector<int>(num_cities));

        // create a vector to store the distances for each path
        vector<double> distances(num_ants);

        // create the ants and let them find a path
        for (int ant = 0; ant < num_ants; ant++)
        {
            // create a vector to store the visited cities for the current ant
            vector<int> visited_cities;

            // choose a random starting city
            uniform_int_distribution<int> dist(0, num_cities - 1);
            int current_city = dist(rng);
            visited_cities.push_back(current_city);

            // move the ant to the next city
            for (int step = 1; step < num_cities; step++)
            {
                // choose the next city using the probability distribution
                pair<int, double> next_city = choose_next_city(current_city, visited_cities, pheromones, cities);
                int city = next_city.first;
                double prob = next_city.second;

                // update the pheromones on the edge between the current city and the chosen city
                pheromones[current_city][city] *= (1.0 - evap_rate);
                pheromones[city][current_city] *= (1.0 - evap_rate);

                double contribution = 1.0 /


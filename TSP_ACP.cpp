#include <fstream>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>
#include <chrono>
#include <ctime>
#include <cmath>
#include <algorithm>
using namespace std;
using namespace std::chrono;

// Intial Values for the Ant Colony (Globa Values)
int numberOfAnts = 25;        // Number of ants                 15
int numberOfIterations = 400; // Number of iterations          1000
double alpha = 0.8;           // Importance of pheromone        1
double beta = 0.7;            // Importance of heuristic information    1
double rho = 0.95;            // Evaporation rate                         1
double Q = 0.3;               // Constant                               1
double tau0 = 0.1;            // Initial pheromone value              0.8

struct problemInstance
{
    string NAME;
    string TYPE;
    string COMMENT;
    int dimension;
    string EDGE_WEIGHT_TYPE;
    int BEST_KNOWN;
    vector<vector<double>> ProblemVector;
    vector<vector<int>> distanceMatrix;
    vector<int> allCities;
    vector<vector<double>> pheromoneMatrix;
};

problemInstance readingTheData(string file)
{
    struct problemInstance instance;
    vector<vector<double>> ProblemVector;
    ifstream infile;
    infile.open(file);
    int linePossition = 0;
    string x;
    string axis;
    if (infile.is_open())
    {
        while (getline(infile, x))
        {
            if (x.find("NAME:") != string::npos)
            {
                int post = x.find(":");
                instance.NAME = x.substr(post + 1);
                linePossition++;
                // remove the end line character
            }
            else if (x.find("TYPE:") != string::npos)
            {
                int post = x.find(":");
                instance.TYPE = x.substr(post + 1);
                linePossition++;
            }
            else if (x.find("COMMENT:") != string::npos)
            {
                int post = x.find(":");
                instance.COMMENT = x.substr(post + 1);
                linePossition++;
            }
            else if (x.find("DIMENSION:") != string::npos)
            {
                int post = x.find(":");
                instance.dimension = stoi(x.substr(post + 1));
                linePossition++;
            }
            else if (x.find("EDGE_WEIGHT_TYPE:") != string::npos)
            {
                int post = x.find(":");
                instance.EDGE_WEIGHT_TYPE = x.substr(post + 1);
                linePossition++;
            }
            else if (x.find("BEST_KNOWN:") != string::npos)
            {
                int post = x.find(":");
                instance.BEST_KNOWN = stoi(x.substr(post + 1));
                linePossition++;
            }
            else if (x.find("NODE_COORD_SECTION") != string::npos)
            {
                linePossition++;
            }
            else if (linePossition > 6)
            {
                vector<double> Cities;
                int i = 0;
                stringstream ss(x);
                string token;
                // Print the Token
                while (ss >> token)
                {
                    if (i != 0 && token != " " && token != "EOF")
                    {
                        Cities.push_back(stod(token));
                    }
                    else
                    {
                        i++;
                    }
                }
                // cout << Cities[0] << " " << Cities[1] << endl;
                ProblemVector.push_back(Cities);
            }
        }
    }
    infile.close();
    instance.ProblemVector = ProblemVector;
    return instance;
}

int enclidianDistanceOfTwoCities(vector<double> city1, vector<double> city2)
{
    double xi = city1[0];
    double yi = city1[1];
    double xj = city2[0];
    double yj = city2[1];
    double distance = sqrt(pow((xi - xj), 2) + pow((yi - yj), 2));
    int roundedValues = round(distance);
    return roundedValues;
}

problemInstance creatingDistanceVector(problemInstance instance)
{
    vector<vector<int>> distanceMatrix(instance.ProblemVector.size(), vector<int>(instance.ProblemVector.size()));
    // cout << "size of the distance matrix: " << distanceMatrix.size() << endl;
    // cout << "size of the Problem Vector: " << instance.dimension << endl;
    for (int i = 0; i < instance.dimension; i++)
    {
        for (int j = 0; j < instance.dimension; j++)
        {
            // cout << j << endl;
            if (i == j)
            {
                // cout << "i == j" << endl;
                distanceMatrix[i][j] = INFINITY;
                // cout << "i == j" << endl;
            }
            else
            {
                distanceMatrix[i][j] = enclidianDistanceOfTwoCities(instance.ProblemVector[i], instance.ProblemVector[j]);
            }
        }
    }
    instance.distanceMatrix = distanceMatrix;
    return instance;
}

// Finding the Cost of the Solution
int costOfTheSolution(vector<int> solution, problemInstance instance)
{
    int cost = 0;
    for (int i = 0; i < solution.size() - 1; i++)
    {
        cost += instance.distanceMatrix[solution[i]][solution[i + 1]];
    }
    return cost;
}

vector<int> allCities(int numberOfCities)
{
    vector<int> cities;
    for (int i = 0; i < numberOfCities; i++)
    {
        cities.push_back(i);
    }
    return cities;
}

vector<int> findingTheNextCity(vector<int> cities, problemInstance instance)
{
    vector<int> nextCities;

    // Fininding all the next Cities
    for (int i = 0; i < instance.allCities.size(); i++)
    {
        for (int j = 0; j < cities.size(); j++)
        {
            if (instance.allCities[i] == cities[j])
            {
                break;
            }
            else if (j == cities.size() - 1)
            {
                nextCities.push_back(instance.allCities[i]);
            }
        }
    }
    return nextCities;
}

vector<double> findingTheProbabilityOfTheUnvisitedCity(vector<int> cities, problemInstance instance, int currentCity)
{
    vector<double> probabilityOfTheUnvisitedCity;

    // Numinator of the probability
    for (int i = 0; i < cities.size(); i++)
    {
        double probability = pow(instance.pheromoneMatrix[currentCity][cities[i]], alpha) * pow((1.0 / instance.distanceMatrix[currentCity][cities[i]]), 5);
        probabilityOfTheUnvisitedCity.push_back(probability);
    }

    // Denominator of the probability
    double denominator = 0;
    for (int i = 0; i < probabilityOfTheUnvisitedCity.size(); i++)
    {
        denominator += probabilityOfTheUnvisitedCity[i];
    }

    // Computing the Probability
    for (int i = 0; i < probabilityOfTheUnvisitedCity.size(); i++)
    {
        probabilityOfTheUnvisitedCity[i] = probabilityOfTheUnvisitedCity[i] / denominator;
    }
    return probabilityOfTheUnvisitedCity;
}

vector<int> computingTheRouletteWheel(vector<double> probabilityOfUnvisitedCities, vector<int> nextCities, vector<int> solution)
{
    // Computing a randome value from 0 to 1
    double randomValue = (double)rand() / RAND_MAX;
    double sum = 0;
    int selectedCity = 0;
    int index = 0;

    for (int i = 0; i < nextCities.size(); i++)
    {
        sum += probabilityOfUnvisitedCities[i];
        if (sum >= randomValue)
        {

            // cout << "Entering Here" << index << endl;
            selectedCity = nextCities[i];
            index = 0;
            break;
        }
        else
        {

            index++;
            // cout << "Not Entering Here" << endl;
        }
    }
    solution.push_back(selectedCity);
    return solution;
}

// Implementing 2-opt for the solution
vector<int> twoOpt(vector<int> solution, problemInstance instance, int numberOfIterations)
{
    // if it goes above 70 seconds then return the solution
    clock_t start = clock();
    int costOfTheSolutionBefore2Opt = costOfTheSolution(solution, instance);
    int costOfTheSolutionAfter2Opt = 0;
    int numberOfIterationsCurrent = 0;
    while (costOfTheSolutionBefore2Opt != costOfTheSolutionAfter2Opt)
    {
        costOfTheSolutionBefore2Opt = costOfTheSolution(solution, instance);
        for (int i = 0; i < solution.size() - 1; i++)
        {
            numberOfIterationsCurrent++;
            for (int j = i + 1; j < solution.size(); j++)
            {
                vector<int> newSolution = solution;
                reverse(newSolution.begin() + i, newSolution.begin() + j);
                costOfTheSolutionAfter2Opt = costOfTheSolution(newSolution, instance);
                if ( costOfTheSolutionAfter2Opt < 0) {
                    return solution;
                }
                if (costOfTheSolutionAfter2Opt < costOfTheSolutionBefore2Opt && costOfTheSolutionAfter2Opt > 0 && costOfTheSolutionAfter2Opt >= instance.BEST_KNOWN )
                {
                    solution = newSolution;
                    costOfTheSolutionBefore2Opt = costOfTheSolutionAfter2Opt;
                    cout << "Cost of the Solution After 2-opt: " << costOfTheSolutionAfter2Opt << " Iteration: " << numberOfIterationsCurrent << " Iterations Given: " << numberOfIterations << endl;
                }
                
            }
            auto end = clock();
            double elapsed_secs = double(end - start) / CLOCKS_PER_SEC;
            if (elapsed_secs > 70)
            {
                return solution;
            }
            if (numberOfIterationsCurrent > numberOfIterations)
            {
                return solution;
            }
        }
        numberOfIterationsCurrent++;
    }
    cout << "Number of Iterations for 2-opt: " << numberOfIterationsCurrent << endl;
    return solution;
}

double AntColonyImplementation(problemInstance problemInstance, vector<int> neightborhoodSolution, int numberOfIterationsForAntColony)
{

    // Creating the vector of all the cities
    vector<int> cities = allCities(problemInstance.dimension);
    problemInstance.allCities = cities;
    int numberOfIterationsForAntColonyT = 0;

    // Creating the pheromone matrix
    vector<vector<double>> pheromoneMatrix(problemInstance.dimension, vector<double>(problemInstance.dimension));
    for (int i = 0; i < problemInstance.dimension; i++)
    {
        
        for (int j = 0; j < problemInstance.dimension; j++)
        {
            if (i == j)
            {
                pheromoneMatrix[i][j] = 0;
            }
            else
            {
                pheromoneMatrix[i][j] = tau0;
            }
        }
    }
    problemInstance.pheromoneMatrix = pheromoneMatrix;

    vector<int> finalSolution;
    double bestCost = costOfTheSolution(neightborhoodSolution, problemInstance);
    vector<int> bestSolution = neightborhoodSolution;
    clock_t start = clock();
    for (int i = 0; i < numberOfIterations; i++)
    {
        // cout << "Iteration: " << i << endl;
        vector<vector<int>> solutions;
        vector<int> costs;
        numberOfIterationsForAntColonyT++;
        

        for (int j = 0; j < numberOfAnts; j++)
        {
            // Choosing a randomg city as the starting point
            int startingCity = rand() % problemInstance.dimension;

            // Creating the solution vector with its path
            vector<int> solution;
            solution.push_back(problemInstance.allCities[startingCity]);
            while (solution.size() != problemInstance.dimension)
            {

                int currentCity = solution.back();

                // Finding all the possible (not visited) cities to go to
                vector<int> nextCities = findingTheNextCity(solution, problemInstance);

                // Finding the probability for each city
                vector<double> probabilityOfCities;
                probabilityOfCities = findingTheProbabilityOfTheUnvisitedCity(nextCities, problemInstance, currentCity);

                // Computing the Roulette Wheel
                solution = computingTheRouletteWheel(probabilityOfCities, nextCities, solution);
            }
            // Adding the first city to the end of the solution to come back to the starting point
            solution.push_back(problemInstance.allCities[startingCity]);

            // Computing the cost of the solution
            int solutionCost = costOfTheSolution(solution, problemInstance);

            // Adding the solution and its cost to the solutions and costs vectors
            solutions.push_back(solution);
            costs.push_back(solutionCost);

            if (solutionCost < bestCost)
            {

                bestCost = solutionCost;
                bestSolution = solution;
                // Print the Best Cost
                cout << "Best Cost: " << bestCost << " for iteration:" << i << endl;
                numberOfIterationsForAntColonyT = 0;
            }
            // cout << "Best Cost: " << bestCost << " for iteration:" << i << endl;

            finalSolution.push_back(bestCost);
        }
        //bestSolution = twoOpt(bestSolution, problemInstance);
        bestCost = costOfTheSolution(bestSolution, problemInstance);

        // Print out the Solution:
        // cout << "Best Solution: " << endl;
        // // Print it in a text file
        // ofstream myfile;
        // myfile.open("bestSolution.txt");
        // for (int i = 0; i < bestSolution.size(); i++)
        // {
        //     myfile << bestSolution[i] << " ";
        // }
        // cout << endl;
        // Times the pheromone by the evaporation rate
        for (int i = 0; i < problemInstance.dimension; i++)
        {
            for (int j = 0; j < problemInstance.dimension; j++)
            {
                pheromoneMatrix[i][j] = pheromoneMatrix[i][j] * rho;
            }
        }
        int length = min(solutions.size(), costs.size());
        for (int i = 0; i < length; i++)
        {
            vector<int> solution = solutions[i];
            double cost = costs[i];
            for (int j = 0; j < problemInstance.dimension; j++)
            {
                // cout << Q/ cost << endl;
                pheromoneMatrix[solution[j]][solution[j + 1]] += Q / cost;
            }
        }

        for (int i = 0; i < problemInstance.dimension; i++)
        {
            pheromoneMatrix[bestSolution[i]][bestSolution[i + 1]] += Q / bestCost;
        }
        problemInstance.pheromoneMatrix = pheromoneMatrix;

        if (numberOfIterationsForAntColonyT == 150)
        {
            break;
        }
        auto stop = high_resolution_clock::now();
        auto duration = (double)(clock() - start) / CLOCKS_PER_SEC;
        //cout << "Time: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
        if (duration > 100)
        {
            break;
        }
    }

    // Print the Pheromone Matrix in a txt file

    problemInstance.pheromoneMatrix = pheromoneMatrix;
    ofstream myfile;
    myfile.open("pheromoneMatrix.txt");
    for (int i = 0; i < problemInstance.dimension; i++)
    {
        for (int j = 0; j < problemInstance.dimension; j++)
        {
            myfile << problemInstance.pheromoneMatrix[i][j] << " ";
        }
        myfile << endl;
    }

    bestSolution = twoOpt(bestSolution, problemInstance, numberOfIterationsForAntColony);
    bestCost = costOfTheSolution(bestSolution, problemInstance);

    // Print the Best Cost
    cout << "Best Cost: " << bestCost << endl;

    // Computing the Gap
    double gap = ((bestCost - problemInstance.BEST_KNOWN) / problemInstance.BEST_KNOWN) * 100;
    cout << "Gap: " << gap << endl;
    return 0;
}

vector<int> nearestNeighbor(problemInstance problemInstance)
{

    vector<int> solution;
    int startingCity = rand() % problemInstance.dimension;

    vector<int> cities = allCities(problemInstance.dimension);
    problemInstance.allCities = cities;

    cout << "Problem Dimension: " << problemInstance.dimension << endl;
    cout << "All Cities Dimension: " << problemInstance.allCities.size() << endl;

    solution.push_back(problemInstance.allCities[startingCity]);
    cout << "Debugging" << endl;
    while (solution.size() != problemInstance.dimension)
    {
        int currentCity = solution.back();
        vector<int> nextCities = findingTheNextCity(solution, problemInstance);
        int nextCity = nextCities[0];
        solution.push_back(nextCity);
    }
    solution.push_back(problemInstance.allCities[startingCity]);
    return solution;
}
// Print it out in a txt file
void printSolution(vector<int> solution, string fileName, problemInstance problemInstance, double time)
{
    // remove the text type from the file name
    // cout << "Printing the solution" << endl;
    string finalFilename = fileName.substr(0, fileName.size() - 4);
    finalFilename = finalFilename.substr(8, finalFilename.size() - 8);
    ofstream myfile;
    // cout << "/results/" + finalFilename + ".txt" << endl;
    myfile.open("/results/" + finalFilename + ".txt");
    myfile << "NAME : " << fileName << endl;
    int cost = costOfTheSolution(solution, problemInstance);
    int gap = ((cost - problemInstance.BEST_KNOWN) / problemInstance.BEST_KNOWN) * 100;
    myfile << "COST: " << cost << endl;
    myfile << "GAP: " << gap << endl;
    myfile << "TIME: " << time << endl;
}

int main()
{
    int seed = 69; // For the memes
    srand(seed);
    vector<string> filesName = {"./data/ch130.tsp", "./data/d198.tsp", "./data/eil76.tsp", "./data/fl1577.tsp", "./data/kroA100.tsp", "./data/lin318.tsp", "./data/pcb442.tsp", "./data/pr439.tsp", "./data/rat783.tsp", "./data/u1060.tsp"};
    uniform_real_distribution<> dis(0, 1);
    mt19937 gen(100);
    default_random_engine re;
    int iterations = 0;
    for (string file : filesName)
    {
        problemInstance instance = readingTheData(file);
        // Print the problem instance
        cout << "NAME: " << instance.NAME << endl;
        cout << "TYPE: " << instance.TYPE << endl;
        cout << "COMMENT: " << instance.COMMENT << endl;
        cout << "DIMENSION: " << instance.dimension << endl;
        cout << "EDGE_WEIGHT_TYPE: " << instance.EDGE_WEIGHT_TYPE << endl;
        cout << "BEST_KNOWN: " << instance.BEST_KNOWN << endl;
        // Print the problem vector
        // for (int i = 0; i < instance.ProblemVector.size(); i++)
        // {
        //     for (int j = 0; j < instance.ProblemVector[i].size(); j++)
        //     {
        //         cout << instance.ProblemVector[i][j] << " ";
        //     }
        //     cout << endl;
        // }
        // return 0;
      
        cout << "Number of Ants: " << numberOfAnts << endl;
        cout << "Number of Iterations: " << numberOfIterations << endl;
        vector<vector<int>> vect{{2, 3}, {4, 5}, {6, 7}};
        instance = creatingDistanceVector(instance);

        // Testing out the Distance Matrix
        vector<vector<double>> tempPoints = {{2, 3}, {4, 5}, {6, 7}};

        struct problemInstance tempInstance;
        tempInstance.ProblemVector = tempPoints;
        tempInstance.dimension = 3;
        tempInstance = creatingDistanceVector(tempInstance);
        // Print the distance matrix to txt file
        ofstream myfile;
        myfile.open("distanceMatrix.txt");
        for (int i = 0; i < instance.dimension; i++)
        {
            for (int j = 0; j < instance.dimension; j++)
            {
                myfile << instance.distanceMatrix[i][j] << " ";
            }
            myfile << endl;
        }

        clock_t start = clock();
        int numberOfIterationsForAntColony = 500;
        // Call the Nearest Neighbor Implementation
        vector<int> currentSolution = nearestNeighbor(instance);

          if (iterations == 3 || iterations == 9 || iterations == 8)
        {
            cout << "Running the Nearest Neighbor Algorithm" << endl;
            currentSolution = twoOpt(currentSolution, instance, 50);
            numberOfAnts = 5;
            numberOfIterations = 5;
            numberOfIterationsForAntColony = 50;
        } 
        else
        {
            numberOfAnts = 25;
            numberOfIterations = 400;
            numberOfIterationsForAntColony = 500;
            
        }

        // vector<int> bestSolution = twoOpt(currentSolution, instance);
        // double bestCost = costOfTheSolution(bestSolution, instance);
        // cout << "Best Cost: " << bestCost << endl;

        // Call the Ant Colony Implementation
        AntColonyImplementation(instance, currentSolution, numberOfIterationsForAntColony);

        auto stop = high_resolution_clock::now();
        auto duration = (double)(clock() - start) / CLOCKS_PER_SEC;
        //cout << "Time: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
        cout << endl;
        printSolution(currentSolution, file, instance, duration);
        // Give the current time in minutes
        // Print the distance matrix
        //  for (int i = 0; i < instance.distanceMatrix.size(); i++)
        //  {
        //      cout << endl;
        //      for (int j = 0; j < instance.distanceMatrix[i].size(); j++)
        //      {
        //          cout << instance.distanceMatrix[i][j] << " ";
        //      }

        // }

        // An example to check if distance matrix is correct

        // Print the cost of the solution
        // vector<int> solution{0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68, 69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 0};
        // cout << "Cost of the solution: " << costOfTheSolution(solution, instance) << endl;
        iterations++;
    }

    return 0;
}
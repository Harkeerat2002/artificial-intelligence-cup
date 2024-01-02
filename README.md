# Traveling Salesman Problem with Ant Colony Optimization

This project implements a solution to the Traveling Salesman Problem (TSP) using Ant Colony Optimization (ACO). The project also uses the Nearest Neighbor heuristic as an initial starting point, and applies a 2-opt local search algorithm to improve the solution.

## Implementation Details

The project is implemented in C++ and consists of the following main components:

- **Nearest Neighbor heuristic**: This is used to get an initial feasible solution to the TSP.
- **Ant Colony Optimization**: This is the main algorithm used to improve upon the initial solution provided by the Nearest Neighbor heuristic.
- **2-opt local search**: This is used to further refine the solution obtained from the ACO.

The project also includes a mechanism to handle larger problems by applying the 2-opt algorithm right after the Nearest Neighbor heuristic. This provides a better starting point for the ACO.

## Running the Code

To compile and run the code, use the following command:

```sh
g++ -Ofast ./TSP_ACP.cpp -o tsp_acp && ./tsp_acp
```

This will generate an executable file tsp_acp which when run, will execute the TSP solution.

## Data Files
The project includes several .tsp files in the data directory. These files contain the data for various instances of the TSP.

## Results
The results of the algorithm are written to text files in the results directory. Each file corresponds to a specific instance of the TSP.

## Future Improvements
The current implementation can be improved in several ways:

- Fixing the 2-opt implementation to achieve better results.
- Running the 2-opt algorithm multiple times to get a better starting point.
- Improving the parameters of the ACO to get better results.


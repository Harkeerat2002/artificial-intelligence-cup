#include <iostream>
#include <algorithm>
#include <vector>
#include <random>
#include <cmath>
#include <fstream>

using namespace std;

struct City
{
    int x;
    int y;
};

// Structure for the TSP data
struct TspData
{
    std::string name;
    std::string comment;
    std::string type;
    int dimension;
    std::string edge_weight_type;
    int best_known;
    std::vector<std::pair<int, std::pair<int, int>>> node_coord_section;
};

// parse the input data and extract the necessary information
void parse_input(vector<City> &cities, int &num_cities, double &best_known)
{
    // Read the Input Data and extract necessary Information
    cin >> num_cities;
    for (int i = 0; i < num_cities; i++)
    {
        int x, y;
        cin >> x >> y;
        cities.push_back({x, y});
    }
}

TspData readTspData(const std::string &fileName)
{
    ifstream infile;
    infile.open("eil76.tsp");
    cout << "Reading file: " << fileName << endl;
    if (infile.is_open())
    {
        string x;

        while (infile >> x ) {
            cout << x << endl;
        }

        infile.close();
    }
    return TspData();
}

int main()
{
    TspData data = readTspData("eli76.tsp");
    // cout << "Name: " << data.name << std::endl;
    // cout << "Comment: " << data.comment << std::endl;
    // cout << "Type: " << data.type << std::endl;
    // cout << "Dimension: " << data.dimension << std::endl;
    // cout << "Edge weight type: " << data.edge_weight_type << std::endl;
    // cout << "Best known: " << data.best_known << std::endl;
    // cout << "Node coord section:" << std::endl;
    // for (const auto &node : data.node_coord_section)
    // {
    //     cout << "  Node " << node.first << ": (" << node.second.first << ", " << node.second.second << ")" << endl;
    // }
    return 0;
}

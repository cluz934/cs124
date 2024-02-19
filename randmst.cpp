// ProgSet 1 - Randomized Minimum Spanning Tree CS124
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>

// Define the adjacency list type for the graph
// Outer index is the vertex, inner pair is the edge and its weight
using Graph = std::vector<std::vector<std::pair<int, double>>>;


// Function to build a basic random graph with edges of random weights distributed N(0, 1)
Graph buildBasicRandGraph(int numOfVertices) {
    Graph randGraph(numOfVertices);
    std::random_device rd; // Non-deterministic seed
    std::mt19937 gen(rd()); // Mersenne Twister engine
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate edges with random weights
    for (int i = 0; i < numOfVertices; ++i) {
        for (int j = i + 1; j < numOfVertices; ++j) {
            double weight = dist(gen); 

            // Add edge i->j
            randGraph[i].push_back({j, weight}); 

            // Since undirected, also add j->i
            randGraph[j].push_back({i, weight}); 
        }
    }

    return randGraph;
}

// Use Kruskal's to find the minimum spanning tree
Graph kruskal(const Graph& graph) {

    std::vector<std::pair<double, std::pair<int, int>>> edges;

    // Iterate through the graph and add all edges to the vector
    for (int i = 0; i < graph.size(); ++i) {
        for (const auto& edge : graph[i]) {
            edges.push_back({edge.second, {i, edge.first}});
        }
    }

    // Sort the edges by weight
    std::sort(edges.begin(), edges.end());

    Graph mst(graph.size());

    // Create a vector to store the parent of each vertex
    std::vector<int> parent(graph.size());

    // Initialize the parent vector
    for (int i = 0; i < graph.size(); ++i) {
        parent[i] = i;
    }

    // Iterate through the edges
    for (const auto& edge : edges) {
        double weight = edge.first;
        int u = edge.second.first;
        int v = edge.second.second;

        // Find the parent of u and v
        int pu = u;
        while (pu != parent[pu]) {
            pu = parent[pu];
        }
        int pv = v;
        while (pv != parent[pv]) {
            pv = parent[pv];
        }

        // If the parents are not the same, add the edge to the minimum spanning tree
        if (pu != pv) {
            mst[u].push_back({v, weight});
            mst[v].push_back({u, weight});
            parent[pu] = pv;
        }
    }

    return mst;
}


// Get the total weight of the minimum spanning tree; takes in the MST of graph G
double getMSTWeight(const Graph& mst) {
    double weight = 0;
    for (int i = 0; i < mst.size(); ++i) {
        for (const auto& edge : mst[i]) {
            weight += edge.second;
        }
    }
    return weight / 2;
}

// Usage: ./randmst 0 numpoints numtrials dimension
int main(int argc, char* argv[]) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " 0 numpoints numtrials dimension" << std::endl;
        return 1;
    }

    int numpoints = std::stoi(argv[2]);
    int numtrials = std::stoi(argv[3]);
    int dimension = std::stoi(argv[4]);

    double totalWeight = 0;

    for (int i = 0; i < numtrials; ++i) {
        Graph randGraph = buildBasicRandGraph(numpoints);
        Graph mst = kruskal(randGraph);
        double mstWeight = getMSTWeight(mst);
        std::cout << "Trial: " << i+1 << std::endl;
        std::cout << "weight = " << mstWeight << std::endl;
        std::cout << std::endl;
        totalWeight += mstWeight;
    }

    //std::cout << numpoints << " " << numtrials << " " << dimension << " " << totalWeight / numtrials << std::endl;

    return 0;
}
// ProgSet 1 - Randomized Minimum Spanning Tree CS124
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>


double k(int n);

// Non-deterministic seed
std::random_device rd; 

// Mersenne Twister engine
std::mt19937 gen(rd()); 

// Define the adjacency list type for the graph
// Outer index is the vertex, inner pair is the edge and its weight
using Graph = std::vector<std::vector<std::pair<int, double>>>;

// Define the Point type for the coordinates of the vertices
using Point = std::vector<double>;

Point generateRandomPoint(std::mt19937& gen, int dimension) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    Point coordinates(dimension);
    for (int i = 0; i < dimension; ++i) {
        coordinates[i] = dist(gen);
    }
    return coordinates;
}

// Function to calculate the Euclidean distance between two points
double EucDist(const Point& p1, const Point& p2) {
    // Ensure the points are of the same dimension
    assert(p1.size() == p2.size());
    double sum = 0;
    for (int i = 0; i < p1.size(); ++i) {
        sum += (p1[i] - p2[i]) * (p1[i] - p2[i]);
    }
    return std::sqrt(sum);
}


// Graph Type 2 and 3
// Function to build a random graph with vertices in N dimension 
// ex: 3D corresponds to unit cube, 4D corresponds to unit hypercube, etc. 
Graph buildRandGraph(int numOfVertices, int dimension) {

    // Initialize the output graph
    Graph randGraph(numOfVertices);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate random points
    std::vector<Point> points(numOfVertices);
    for (int i = 0; i < numOfVertices; ++i) {
        points[i] = generateRandomPoint(gen, dimension);
    }

    // Generate edges with weights equal to the Euclidean distance between the points
    for (int i = 0; i < numOfVertices; ++i) {
        for (int j = i + 1; j < numOfVertices; ++j) {
            double weight = EucDist(points[i], points[j]);

            // Add edge i->j
            randGraph[i].push_back(std::make_pair(j, weight)); 

            // Since undirected, also add j->i
            randGraph[j].push_back(std::make_pair(i, weight)); 
        }
    }

    return randGraph;
}


// Graph Type 1
// Function to build a basic random graph with edges of random weights distributed N(0, 1)
// For memory efficiency, we prune as we go, so we only store the edges with weights less than k(n)
Graph buildBasicRandGraph(int numOfVertices) {

    // Initialize the output graph
    Graph randGraph(numOfVertices);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate edges with weights distributed N(0, 1)
    for (int i = 0; i < numOfVertices; ++i) {
        for (int j = i + 1; j < numOfVertices; ++j) {
            double weight = dist(gen);

            // Prune the edge if the weight is greater than k(n)
            // Prune as we go to save memory
            if (weight < k(numOfVertices)) {
                // Add edge i->j
                randGraph[i].push_back(std::make_pair(j, weight)); 

                // Since undirected, also add j->i
                randGraph[j].push_back(std::make_pair(i, weight));
            }
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
            edges.push_back(std::make_pair(edge.second, std::make_pair(i, edge.first)));
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
            mst[u].push_back(std::make_pair(v, weight));
            mst[v].push_back(std::make_pair(u, weight));
            parent[pu] = pv;
        }
    }

    return mst;
}


double k(int n) {
    // TODO: Find a better function for k(n)
    return 3*std::log(n) / n;
}

// PRUNING THE TREES
// To handle large number of vertices, we can prune a tree before running Kruskal's
// The minimum spanning tree is extremely unlikely to use any edge of weight greater than k(n) for some function k(n)
// We can estimate k(n) using small values of n, and then try to throw away edges of weight larger than k(n) as we increase the input size
// throwing away too many edges may cause problems, but throwing away too few edges may not help at all
Graph pruneTree(const Graph& graph, int n) {
    // Output graph
    Graph prunedGraph(graph.size());

    // threshold is the value of k(n), edges with weight greater than threshold are pruned
    double threshold = k(n);

    // Iterate through the graph and add all edges to the pruned graph
    for (int i = 0; i < graph.size(); ++i) {
        for (const auto& edge : graph[i]) {
            if (edge.second <= threshold) {
                prunedGraph[i].push_back(edge);
            }
        }
    }
    return prunedGraph;
}

// Get largest edge weight in the MST of a graph
// Takes in the MST of graph G
double getLargestEdgeWeight(const Graph& mst) {
    double largestWeight = 0;
    for (int i = 0; i < mst.size(); ++i) {
        for (const auto& edge : mst[i]) {
            if (edge.second > largestWeight) {
                largestWeight = edge.second;
            }
        }
    }
    return largestWeight;
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

    double totalLength = 0;

    // Run the trials
    for (int i = 0; i < numtrials; ++i) {
        // Testing Graph Creation in n-dimensional space given by "dimension"
        Graph graph;
        Graph mst;
        if (dimension == 0) {
            graph = buildBasicRandGraph(numpoints);
            mst = kruskal(graph);
        } else {
            graph = buildRandGraph(numpoints, dimension);
            Graph prunedTree = pruneTree(graph, numpoints);
            mst = kruskal(prunedTree);
        }

        double mstWeight = getMSTWeight(mst);
        totalLength += mstWeight;
        std::cout << "Trial " << i+1 << " MST Weight: " << mstWeight << std::endl;
    }
    std::cout << std::endl;
    std::cout << "Number of Vertices: " << numpoints << std::endl;
    std::cout << "Number of Trials: " << numtrials << std::endl;
    std::cout << "Dimension: " << dimension << std::endl;
    std::cout << "Average MST Weight: " << totalLength / numtrials << std::endl;
    return 0;
}
#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include <chrono>
#include <future>
#include <thread>

// Non-deterministic seed
std::random_device rd; 

// Mersenne Twister engine
std::mt19937 gen(rd()); 

// Threshold function
double k(int n, int dimension);

using Point = std::vector<double>;

// Edge struct for the graph
// Handles both 0D and higher dimensions


using Graph = std::vector<std::vector<std::pair<int, double>>>;


// We need union-find to find the minimum spanning tree via Kruskal's
class UnionFind {
    private:
        std::vector<int> parent;
        std::vector<int> rank;

    public:

    // Constructor
    UnionFind(int N): parent(N), rank(N, 0) {
        for (int i = 0; i < N; ++i) {
            parent[i] = i;
        }
    }

    int find(int u) {
        if (u != parent[u]) {
            // Path compression
            parent[u] = find(parent[u]);
        }
        return parent[u];
    }

    // Unite the sets that contain u and v
    void unite(int u, int v) {
        int pu = find(u);
        int pv = find(v);

        // If the parents are the same, then the vertices are already in the same set
        if (pu == pv) {
            return;
        }

        // Union by rank
        // If the rank of pu is less than the rank of pv, then make pv the parent of pu
        if (rank[pu] < rank[pv]) { 
            parent[pu] = pv;
        } else if (rank[pu] > rank[pv]) { // Rank of pv is less than rank of pu
            parent[pv] = pu;
        } else { // Rank of pv is equal to rank of pu
            parent[pu] = pv;
            rank[pv]++;
        }
    }
};

// Threshold function
double k(int n, int dimension) {

    // TODO: Find a better function for k(n) in higher dimensions
    if (dimension == 0) {
        return 3*std::log2(n) / n;
    } else if (dimension == 2) {
        return 1;
    } else if (dimension == 3) {
        return 1;
    } else if (dimension == 4) {
        return 1;
    } else {
        return 1;
    }
}

// Generate a random point in n-dimensional space
Point generateRandomPoint(int dimension, std::mt19937& gen) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    Point point(dimension);
    for (int i = 0; i < dimension; ++i) {
        point[i] = dist(gen);
    }
    return point;
}

// Euclidean distance between two points in n-dimensional space
// Input: Two points and the dimension
// Output: The Euclidean distance between the two points
double EucDistance(const Point& p1, const Point& p2, int dimension) {
        double sum = 0.0;
        for (int i = 0; i < dimension; ++i) {
            double diff = p1[i] - p2[i];
            sum += diff * diff;
        }
        return std::sqrt(sum);
}


// Graph Type 1: Random Graph in 0 dimensions
// For memory efficency, we prune the edges that have a weight greater than k(n) as we build the graph
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
            if (weight < k(numOfVertices, 0)) {
                // Add edge i->j
                randGraph[i].push_back(std::make_pair(j, weight)); 
            }
        }
    }
    return randGraph;
}


// Graph Type 2, 3: Random Graph in n dimensions (up to 4 dimensions)
// Build a random graph in n-dimensional space; Prune edges that have a weight greater than k(n)
// Output: A vector of edges
Graph buildRand_NDGraph(int numOfVertices, int dimension) {

    // Initialize the output graph
    Graph randGraph(numOfVertices);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    // Generate random points
    std::vector<Point> points(numOfVertices);
    for (int i = 0; i < numOfVertices; ++i) {
        points[i] = generateRandomPoint(dimension, gen);
    }

    // Generate edges with weights equal to the Euclidean distance between the points
    for (int i = 0; i < numOfVertices; ++i) {
        for (int j = i + 1; j < numOfVertices; ++j) {
            double weight = EucDistance(points[i], points[j], dimension);

            if (weight < k(numOfVertices, dimension)) {
                randGraph[i].push_back(std::make_pair(j, weight));
            }
        }
    }

    return randGraph;
}


// Kruskal's algorithm to find the minimum spanning tree using our UnionFind data structure and our edge struct
Graph kruskals(int numOfVertices, const Graph& graph) {
    // Initialize the output graph
    Graph mst(numOfVertices);

    // Initialize the UnionFind data structure
    UnionFind uf(numOfVertices);

    // Sort the edges by weight
    std::vector<std::pair<int, std::pair<int, double>>> edges;
    for (int i = 0; i < numOfVertices; ++i) {
        for (const auto& edge : graph[i]) {
            edges.push_back(std::make_pair(i, edge));
        }
    }
    std::sort(edges.begin(), edges.end(), [](const auto& a, const auto& b) {
        return a.second.second < b.second.second;
    });

    // Run Kruskal's algorithm
    for (const auto& edge : edges) {
        int u = edge.first;
        int v = edge.second.first;
        double weight = edge.second.second;

        if (uf.find(u) != uf.find(v)) {
            mst[u].push_back(std::make_pair(v, weight));
            mst[v].push_back(std::make_pair(u, weight));
            uf.unite(u, v);
        }
    }

    return mst;
}

double getMSTWeight(const Graph& mst) {
    double weight = 0;
    for (int i = 0; i < mst.size(); ++i) {
        for (const auto& edge : mst[i]) {
            weight += edge.second;
        }
    }

    // We divide by 2 because the MST is undirected
    return weight / 2;
}

// Finds longest edge in MST
double getLongestEdge(const Graph& mst) {
    double longest = 0;
    for (int i = 0; i < mst.size(); ++i) {
        for (const auto& edge : mst[i]) {
            if (edge.second > longest) {
                longest = edge.second;
            }
        }
    }
    return longest;
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
    double longestEdge = 0;
    // Generate the edges for the graph and run the trials
    for (int i = 0; i < numtrials; ++i) {
        Graph trial_graph;

        auto start_graph = std::chrono::high_resolution_clock::now();

    
        if (dimension == 0) {
            trial_graph = buildBasicRandGraph(numpoints);
        } else {
            // Replace buildRand_NDGraph with buildRand_NDGraph_MT for multithreaded execution
            trial_graph = buildRand_NDGraph(numpoints, dimension);
        }

        auto end_graph = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed_graph = end_graph - start_graph;

        auto start_kruskal = std::chrono::high_resolution_clock::now();
        Graph mst = kruskals(numpoints, trial_graph);
        auto end_kruskal = std::chrono::high_resolution_clock::now();

        double weight = getMSTWeight(mst);
        double longest = getLongestEdge(mst);

        longestEdge += longest;
        totalLength += weight;

        // Trial Overview
        // std::cout << "Trial: " << i + 1 << std::endl;
        // std::cout << "MST Weight: " << weight << std::endl;
        // std::cout << "Longest Edge: " << longest << std::endl;
        // std::cout << "Time to generate graph: " << elapsed_graph.count() << "s" << std::endl;
        // std::cout << "Time to run Kruskal's: " << std::chrono::duration<double>(end_kruskal - start_kruskal).count() << "s" << std::endl;
        // std::cout << std::endl;
    }

        // std::cout << "Number of Trials: " << numtrials << std::endl;
        // std::cout << "Number of Vertices: " << numpoints << std::endl;
        // std::cout << "Dimension: " << dimension << std::endl;
        // std::cout << "Average MST Weight: " << totalLength / numtrials << std::endl;
        // std::cout << "Average Longest Edge: " << longestEdge / numtrials << std::endl;
    
    // Output format: average numpoints numtrials dimension
    // where average is the average minimum spanning tree weight over the trials.

    std::cout << totalLength / numtrials << " " << numpoints << " " << numtrials << " " << dimension << std::endl;

    return 0;
}
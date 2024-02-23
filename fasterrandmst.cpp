#include <random>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>


// Non-deterministic seed
std::random_device rd; 

// Mersenne Twister engine
std::mt19937 gen(rd()); 

// Threshold function
double k(int n);

struct Edge {
    int u; // Vertex u
    int v; // Vertex v
    double weight; // Weight of the edge

    // Constructor
    Edge(int _u, int _v, double _weight) : u(_u), v(_v), weight(_weight) {}

    // Comparator for sorting
    bool operator<(const Edge& e) const {
        return weight < e.weight;
    }
};

std::vector<Edge> graphEdges;


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

double k(int n) {
    // TODO: Find a better function for k(n)
    return 3*std::log(n) / n;
}

// For memory efficency, we prune the edges that have a weight greater than k(n) as we build the graph
std::vector<Edge> buildBasicRandGraph(int numOfVertices) {
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    std::vector<Edge> edges;

    // Generate all possible edges with random weights
    for (int i = 0; i < numOfVertices; ++i) {
        for (int j = i + 1; j < numOfVertices; ++j) {

            // Generate the random weight
            double weight = dist(gen); 

            // Prune the edge if the weight is greater than k(n)
            if (weight <= k(numOfVertices)) {

                // Add the edge i->j
                edges.emplace_back(i, j, weight);

                // Add the edge j->i
                edges.emplace_back(j, i, weight);

            }
        }
    }
    return edges;
}

// Kruskal's using union-find
// Returns: The minimum spanning tree of graph G
std::vector<Edge> kruskals(int numOfVertices, std::vector<Edge>& edges) {
    std::vector<Edge> mst;
    std::sort(edges.begin(), edges.end());

    UnionFind uf(numOfVertices);

    for (const auto& edge : edges) {
        int u = edge.u;
        int v = edge.v;
        double weight = edge.weight;

        if (uf.find(u) != uf.find(v)) {
            mst.push_back(edge);
            uf.unite(u, v);
        }
    }

    return mst;
}

double getMSTWeight(const std::vector<Edge>& mst) {
    double weight = 0;
    for (const auto& edge : mst) {
        weight += edge.weight;
    }
    return weight;
}


// Usage: ./fasterrandmst 0 numpoints numtrials dimension
int main(int argc, char* argv[]) {

    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " 0 numpoints numtrials dimension" << std::endl;
        return 1;
    }

    int numpoints = std::stoi(argv[2]);
    int numtrials = std::stoi(argv[3]);
    int dimension = std::stoi(argv[4]);

    // Generate the edges for the graph
    graphEdges = buildBasicRandGraph(numpoints);
    double totalLength = 0;

    // Run the trials
    for (int i = 0; i < numtrials; ++i) {
        // Testing Graph Creation in n-dimensional space given by "dimension"
        std::vector<Edge> mst = kruskals(numpoints, graphEdges);
        double weight = getMSTWeight(mst);
        totalLength += weight;
        std::cout << "Trial " << i << ": " << weight << std::endl;
    }
    std::cout << "Average MST Weight: " << totalLength / numtrials << std::endl;

    return 0;
}
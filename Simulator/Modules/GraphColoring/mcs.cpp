#include "mcs.h"
#include <queue>
#include <list>
#include <chrono>       // std::chrono::system_clock
#include <random>

using std::queue;
using std::cout;
using std::cerr;
using std::endl;

using namespace GAIA::GraphColoring;

vector<int>& GAIA::GraphColoring::Mcs::color()
{
    //std::list<int> temp_graph;
    std::vector<int> temp_graph;
    for (size_t i = 0; i < graph.size(); i++)
    {
        temp_graph.push_back(i);
    }


    vector<int> weight(temp_graph.size());
    queue<int> ordering;

    std::cout << "Initializing.\n";

    // Initially set the weight of each node to 0
    for (int i = 0; i < weight.size(); ++i) {
        weight[i] = 0;
    }

    std::cout << "Working through all the nodes in the graph to update maximum weight.\n";

    // Work through all the nodes in the graph, choosing the node
    // with maximum weight, then add that node to the queue. Increase
    // the weight of the queued nodes neighbors by 1. Continue until
    // every node in the graph has been added to the queue

    int percentage = 0;

    std::vector<int> coloringOrder(this->graph.size());
    for (int iNode = 0; iNode < this->graph.size(); iNode++) {
        coloringOrder[iNode] = iNode;

    }


    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine e(seed);
    std::shuffle(std::begin(coloringOrder), std::end(coloringOrder), e);

    for (int i = 0; i < this->graph.size(); i++) {
        int iNode = coloringOrder[i];

        int max_weight = -1;
        int max_vertex = -1;

        // Out of the remaining nodes, find the node with the highest weight
        ;
        int maxWId;
        for (int j = 0; j < temp_graph.size(); j++) {
            int nodeId = temp_graph[j];
            if (nodeId < 0)
            {
                continue;
            }
            if (weight[nodeId] > max_weight) {
                max_weight = weight[nodeId];
                max_vertex = nodeId;
                maxWId = j;
            }
        }
        if (max_vertex == -1) {
            cerr << "Error: Could not find a max weight node in the graph (reason unknown)" << endl;
            this->graph_colors = vector<int>();
            return graph_colors;
        }

        // Add highest weight node to the queue and increment all of its
        // neighbors weights by 1
        ordering.push(max_vertex);
        for (unsigned j = 0; j < graph[max_vertex].size(); j++) {
            weight[graph[max_vertex][j]] += 1;
        }

        // Remove the maximum weight node from the graph so that it won't
        // be accidentally added again
        //temp_graph.erase(graphIterMaxW); // 73010 used 3:21
        //*graphIterMaxW = -1;             // 73010 used 2:31
        temp_graph[maxWId] = -1;                // 73010 used 2:10

        //std::cout << i << "  th iteration.\n";
        if (100.0 * (double)iNode / this->graph.size() > percentage)
        {
            if (verbose)
            {
                std::cout << percentage << "%  finished." << endl;
            }
            percentage += 1;

            // update temp_graph to remove negative ids
            std::vector<int> temp_graph_new;
            for (size_t j = 0; j < temp_graph.size(); j++)
            {
                if (temp_graph[j] > 0)
                {
                    temp_graph_new.push_back(temp_graph[j]);
                }
            }
            temp_graph = std::move(temp_graph_new);

        }

    }

    int sizeOrdering = ordering.size();
    percentage = 0;
    std::cout << "Work through the queue in order and color each node.\n";

    // Work through the queue in order and color each node
    while (!ordering.empty()) {
        int color = 0;

        // Find the lowest possible graph_colors for this node between
        // its neighbors
        int min = ordering.front();

        //Thanks to Michael Kochte @ Universitaet Stuttgart for the below speedup snippit

        //Collect color numbers of neighbors
        vector<int> colorvec;
        for (unsigned i = 0; i < graph[min].size(); i++) {
            int col = graph_colors[graph[min][i]];
            if (std::find(colorvec.cbegin(), colorvec.cend(), col) == colorvec.cend()) {
                colorvec.push_back(col);
            }
        }
        
        //Sort and uniquify
        std::sort(colorvec.begin(), colorvec.end());

        //Pick the lowest color not contained
        int newcolor = 0;
        for (unsigned i = 0; i < colorvec.size(); i++) {
            if (colorvec[i] == newcolor) {
                newcolor++;
            }
        }
        color = newcolor;

        this->graph_colors[min] = color;
        ordering.pop();

        if (100.0 * (double)(sizeOrdering - ordering.size()) / sizeOrdering > percentage)
        {
            if (verbose)
            {
                std::cout << percentage << "%  finished." << endl;;
            }
            percentage += 1;

        }

    }
    return this->graph_colors;
}

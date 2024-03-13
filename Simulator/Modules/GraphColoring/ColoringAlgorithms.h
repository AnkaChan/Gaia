#pragma once

#include "Graph.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>



namespace GAIA {
	namespace GraphColoring {
        using std::map;
        using std::string;
        using std::vector;

        class GraphColor {
        protected:
            /* Internal Members */
            vector<vector<int>> graph; // node x neighbors
            vector<int> graph_colors;
            vector<vector<int>> categories;

            bool verbose = false;
        public:
            typedef std::shared_ptr<GraphColor> SharedPtr;
            typedef GraphColor* Ptr;


            /* Constructors */
            GraphColor();
            GraphColor(const Graph& graph);
            GraphColor(const std::vector<vector<int>>& graph);

            /* Mutators */
            virtual vector<int>& color() = 0;
            void set_graph(const vector<vector<int>>& new_graph) { this->graph = new_graph; }
            void modify_graph(int node, const vector<int>& neighbors) { this->graph[node] = neighbors; }

            /* Accessors */
            virtual string get_algorithm() = 0;
            unsigned size() { return this->graph.size(); }
            bool is_colored();
            vector<int>& get_coloring() { return this->graph_colors; }
            int get_color(int node) { return graph_colors[node]; }
            int get_num_colors();
            bool is_valid();

            /* Print functions */
            void convertToColoredCategories();

            /* Functions for coloring balancing */
            float findLargestSmallestCategories(int& biggestCategory, int& smallestCategory);
            // return the category id of the changable node, not the node id, -1 if not changable 
            int findChangableNodeInCategory(int sourceColor, int destinationColor);
            void changeColor(int sourceColor, int categoryId, int destinationColor);
            bool changable(int node, int destinationColor);
            void balanceColoredCategories(float goalMaxMinRatio = 1.5);
            void saveColoringCategories(std::string outputFile);

        };

        //class GraphClustering : public GraphColor {
        //public:
        //    GraphClustering(const std::vector<vector<int>>& inGraph, const std::vector<vector<float>>& inGraphWeights);
        //    void cluster(int Ks);
        //    // 1 - sharedVerts / (8 - sharedVerts)
        //    vector<vector<float>> graph_weights;
        //    // graph_clusters[i] is the cluster that node i belongs to;
        //    // if graph_clusters[i] returns -1, the constraint has not been placed in a cluster yet;
        //    vector<int> orgGraphNodeClusters;
        //    // the new graph made of the clustered nodes
        //    vector<vector<int>> clusteredGraph;
        //    // the original nodes that each node of the clustered graph contains
        //    vector<vector<int>> clusteredGraphNodes;
        //    int minimumUnclusteredWeightNode();
        //    void saveClusteredColoringCategories(int numColors, vector<int>& clusteredGraphColors, std::string outputFile);

        //private:
        //    virtual vector<int>& color() { std::cout << "This is a clustering graph, don't use it for coloring!\n"; return graph_colors; }
        //    virtual string get_algorithm() { return "clustering"; }

        //};
	}
}
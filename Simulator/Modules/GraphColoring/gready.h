#pragma once
#include "ColoringAlgorithms.h"

namespace GAIA {
	namespace GraphColoring {

		/*
		* Coloring algorithm provided by 
		* Ton-That, Quoc-Minh, Paul G. Kry, and Sheldon Andrews. "Parallel block Neo-Hookean XPBD using graph clustering." Computers & Graphics 110 (2023): 1-10.
		* Faster than MCS but results are inferior
		*/

		class OrderedGreedy : public GraphColor {
		public:
			/* Constructors */
			OrderedGreedy(const Graph& graph) : GraphColor(graph) {}
			OrderedGreedy(vector<vector<int>>& graph) : GraphColor(graph) {}

			virtual int nextNode();
			void reduceDegree(int iNode);
			/* Mutators */
			vector<int>& color();

			/* Accessors */
			string get_algorithm() { return "OrderedGreedy"; }

			std::vector<int> degrees;
			std::vector<bool> colored;
		};
	}
}

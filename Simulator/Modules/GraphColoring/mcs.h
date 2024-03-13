#pragma once

#include "ColoringAlgorithms.h"


namespace GAIA {
	namespace GraphColoring {

		/*
		* MCS (Register Allocation via Coloring of Chordal Graphs - Magno et al.) 
		*/
		class Mcs : public GraphColor {
		public:
			/* Constructors */
			Mcs(const Graph& graph) : GraphColor(graph) {}
			Mcs(vector<vector<int>>& graph) : GraphColor(graph) {}

			/* Mutators */
			vector<int>& color();

			/* Accessors */
			string get_algorithm() { return "MCS"; }
		};
	}
}
#pragma once
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

#include "Graph.h"

namespace GAIA {
	namespace GraphColoring {
		struct TetMeshVertexGraph : Graph
		{
			void fromMesh(void* pMesh) {
				TMeshStaticF::Ptr pTM = (TMeshStaticF::Ptr)pMesh;
				numNodes = pTM->numVertices();
				edges.clear();

				std::set<Edge> edgeSet;

				for (TMeshStaticF::VPtr pV : TIt::TM_VIterator(pTM))
				{
					for (TMeshStaticF::VPtr pNeiV : TIt::V_VIterator(pV))
					{
						Edge e;
						if (pV->id() < pNeiV->id()) {
							e = { pV->id(), pNeiV->id() };
						}
						else
						{
							e = { pNeiV->id(), pV->id() };
						}
						edgeSet.insert(e);
					}
				}

				for (auto& e : edgeSet)
				{
					edges.push_back(e);
				}
			}
		};

	}
}
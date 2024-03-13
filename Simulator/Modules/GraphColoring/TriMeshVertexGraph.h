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
		struct TriMeshVertexGraph : Graph
		{
			TriMeshVertexGraph(bool inAddEdgeBasedBending_in)
				: addEdgeBasedBending(inAddEdgeBasedBending_in) 
			{}
			

			void fromMesh(void* pMesh_in) {
				TriMeshStaticF::Ptr pMesh = (TriMeshStaticF::Ptr)pMesh_in;
				numNodes = pMesh->numVertices();
				edges.clear();

				for (TriMeshStaticF::EPtr pE : It::MEIterator(pMesh))
				{
					edges.push_back({pE->halfedge()->source()->id(), pE->halfedge()->target()->id() });

					if (addEdgeBasedBending)
					{
						TriMeshStaticF::HEPtr pHEDual = pE->halfedge()->he_sym();
						if (pHEDual != nullptr)
						{
							int vIdCross1 = pE->halfedge()->he_next()->target()->id();
							int vIdCross2 = pHEDual->he_next()->target()->id();
							edges.push_back({ vIdCross1, vIdCross2 });
						}
					}
				}

			}
			protected:
			bool addEdgeBasedBending;

		};

	}
}
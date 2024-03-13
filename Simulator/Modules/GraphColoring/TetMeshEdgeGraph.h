#pragma once
#include <vector>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>
#include <map>

namespace GAIA {
	namespace GraphColoring {
		struct TetMeshEdgeGraph : Graph {
			std::vector<float> edgeWeights;

			void fromMesh(void* pMesh) {
				TMeshStaticF::Ptr pTM = (TMeshStaticF::Ptr)pMesh;

				numNodes = pTM->numEdges();
				edges.clear();
				std::set<std::array<int, 2>> edgesSet;

				std::map<TMeshStaticF::EPtr, int> edgeIdMap;

				int eId = 0;
				for (TMeshStaticF::EPtr pE : TIt::TM_EIterator(pTM))
				{
					edgeIdMap.insert(std::pair<TMeshStaticF::EPtr, int>(pE, eId));
					eId++;
				}

				for (TMeshStaticF::VPtr pV : TIt::TM_VIterator(pTM))
				{
					std::vector<int> adjacentEdgeIds;
					for (TMeshStaticF::EPtr pE : TIt::V_EIterator(pV))
					{
						adjacentEdgeIds.push_back(edgeIdMap[pE]);
					}
					std::sort(adjacentEdgeIds.begin(), adjacentEdgeIds.end());
					// all the edge adjacent to this vertex forms a complete graph
					for (size_t iT = 0; iT < adjacentEdgeIds.size(); iT++)
					{
						for (size_t iT2 = iT + 1; iT2 < adjacentEdgeIds.size(); iT2++) {
							std::array<int, 2> edge;
							if (adjacentEdgeIds[iT] >= adjacentEdgeIds[iT2])
								edge = { adjacentEdgeIds[iT], adjacentEdgeIds[iT2] };
							else
								edge = { adjacentEdgeIds[iT2], adjacentEdgeIds[iT] };
							assert(edgesSet.find(edge) == edgesSet.end());

							edgesSet.insert(edge);
							edges.push_back(edge);
							// adjacent edges share exactly 1 common vertex
							edgeWeights.push_back(1.f - 1.f / 3.f);
						}
					}
				}
			};

			void saveColFile(std::string outFile) {
				std::ofstream ofs(outFile);
				if (ofs.is_open())
				{
					ofs << "p edge " << numNodes << " " << edges.size() << "\n";
					for (size_t iEdge = 0; iEdge < edges.size(); iEdge++)
					{
						ofs << "e " << edges[iEdge][0] + 1 << " " << edges[iEdge][1] + 1 << "\n";
					}
					ofs.close();
				}
				else
				{
					std::cout << "Fail to open: " << outFile;
				}

				std::string outWeightsFile = outFile + ".weights.json";

				std::ofstream ofsWeights(outWeightsFile);
				if (ofsWeights.is_open())
				{
					nlohmann::json outWeights = edgeWeights;
					ofsWeights << outWeights.dump();
					ofsWeights.close();
				}
				else
				{
					std::cout << "Fail to open: " << outWeightsFile;
				}
			}
		};
	}
}
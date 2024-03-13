#pragma once
#include <vector>
#include <array>
#include <algorithm>
#include <string>
#include <fstream>
#include <iostream>

#include "Json/json.hpp"

#include <MeshFrame/TetMesh/TMeshStaticLibHeaders.h>


#include "Graph.h"

namespace GAIA {
	namespace GraphColoring {

		struct TetMeshTetGraph : Graph{
			size_t numNodes;

			void fromMesh(void* pMesh) {
				TMeshStaticF::Ptr pTM = (TMeshStaticF::Ptr)pMesh;
				numNodes = pTM->numTets();
				edges.clear();
				std::set<std::array<int, 2>> edgesSet;

				for (TMeshStaticF::VPtr pV : TIt::TM_VIterator(pTM))
				{
					std::vector<int> adjacentTetIds;
					for (TMeshStaticF::TVPtr pTV : TIt::V_TVIterator(pV))
					{
						TMeshStaticF::TPtr pT = TMeshStaticF::TVertexTet(pTV);
						adjacentTetIds.push_back(pT->id());
					}
					std::sort(adjacentTetIds.begin(), adjacentTetIds.end());
					// all the tets adjacent to this vertex forms a complete graph
					for (size_t iT = 0; iT < adjacentTetIds.size(); iT++)
					{
						for (size_t iT2 = iT + 1; iT2 < adjacentTetIds.size(); iT2++) {

							TMeshStaticF::TPtr pT1 = pTM->idTet(adjacentTetIds[iT]);
							TMeshStaticF::TPtr pT2 = pTM->idTet(adjacentTetIds[iT2]);

							int numSharedVerts = 0;
							for (TMeshStaticF::VPtr pVT1 : TIt::T_VIterator(pT1))
							{
								for (TMeshStaticF::VPtr pVT2 : TIt::T_VIterator(pT2)) {
									if (pVT2 == pVT1)
									{
										++numSharedVerts;
										break;
									}
								}
							}
							Edge edge;
							if (adjacentTetIds[iT] <= adjacentTetIds[iT2])
								edge = { adjacentTetIds[iT], adjacentTetIds[iT2] };
							else
								edge = { adjacentTetIds[iT2], adjacentTetIds[iT] };

							if (edgesSet.find(edge) == edgesSet.end()) {
								edges.push_back(edge);
								edgesSet.insert(edge);
								edgeWeights.push_back(1.f - float(numSharedVerts) / (8.f - numSharedVerts));
							}
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
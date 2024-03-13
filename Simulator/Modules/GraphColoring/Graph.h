#pragma once
#include <MeshFrame/TetMesh/TMeshStaticLibHeaders.h>
#include <MeshFrame/TriMesh/MeshStatic.h>

namespace GAIA {
	namespace GraphColoring {
		using TMeshStaticF = MF::TetMesh::TMeshStaticF;
		typedef MF::TetMesh::TIterators<TMeshStaticF> TIt;

		using TriMeshStaticF = MF::TriMesh::TriMeshStaticF;
		typedef MF::TriMesh::TriMeshStaticIteratorF It;

		struct Graph {
			size_t numNodes;
			typedef std::array<int, 2> Edge;
			std::vector<Edge> edges;
			std::vector<float> edgeWeights;

			virtual void fromMesh(void* pMesh) = 0;

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

			}
		};
	}
}


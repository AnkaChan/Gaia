#pragma once
#include "../TetMesh/TetMeshFEM.h"
#include "../TriMesh/TriMesh.h"

namespace GAIA {
	//void writeObj(const char* output, & tIf);
	//void writeAllToObj(const char* output, std::vector<SoftBodyPtr>& softBodies, bool saveAllInOneFile = false);
	void writeAllToPLY(const char* output, std::vector<TetMeshFEM::SharedPtr>& tetMeshes, bool saveAllInOneFile = false, bool addModelName = false);
	void writeAllToBinary(const char* output, std::vector<TetMeshFEM::SharedPtr>& tetMeshes);
	void writeAllToPLY(const char* output, std::vector<TriMeshFEM::SharedPtr>& tetMeshes, bool saveAllInOneFile = false);
	//void writeTetMeshToVtk(const char* output, TetMeshInterFace& tIf);
	//void writeAllTetMeshToVtk(const char* output, std::vector<SoftBodyPtr>& softBodies);
	//void writeAllTetMeshToT(const char* output, std::vector<SoftBodyPtr>& softBodies);
	void saveAsPly(const std::string filePath, const std::vector<std::array<FloatingType, 3>>& verts, const std::vector<std::array<int, 3>>& faces);
	void saveAsPly(const std::string filePath, const std::vector<std::array<FloatingType, 3>>& verts, const std::vector<std::array<int, 3>>& faces);
}

#pragma once

namespace GAIA {

	struct TriMeshEdgeGPU
	{
		int eV1, eV2,
			eV12Next, // the vertex opposite to edge (v1, v2), in the nei triangle whose orientation matches that order
			eV21Next; // the vertex opposite to edge (v2, v1), in the nei triangle whose orientation matches that order

		int fId1,  //  the nei triangle whose orientation matches (v1, v2)
			fId2;  //  the nei triangle whose orientation matches (v2, v1)
	};

	struct TriMeshTopologyGPU final
	{
		int32_t* faceVIds;

		int32_t* vertexNeighborFaces;
		int32_t* vertexNeighborFaces_vertexOrder;
		int32_t* vertexNeighborFaces_infos;

		int32_t* vertexNeighborBendings;
		int32_t* vertexNeighborBendings_vertexOrder;
		int32_t* vertexNeighborBendings_infos;

	};

	struct TriMeshFEMGPU
	{
		typedef TriMeshFEMGPU* Ptr;

		// topology
		TriMeshTopologyGPU* pTopology;

		//* non-variant data
		// data arranged per vertex
		FloatingTypeGPU* vertexMass;

		// data arranged per face
		FloatingTypeGPU* DmInvs;            // face's mats, 4 x nTets
		FloatingTypeGPU* faceRestArea;
		FloatingTypeGPU* faceInvRestArea;
		//* non-variant data

		// flatten vertex data: 3*nVerts
		FloatingTypeGPU* vertPos;
		FloatingTypeGPU* vertPrevPos;
		FloatingTypeGPU* velocity;

		int8_t* vertexFixedMask;

		bool activeForCollision = true;
		bool activeForMaterialSolve = true;
	};

}
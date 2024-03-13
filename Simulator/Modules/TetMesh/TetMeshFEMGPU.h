#pragma once

#define VERTEX_BUFFER_STRIDE 3

#include "CuMatrix/Buffers/ManagedBuffer.h"
typedef float FloatingTypeGPU;
typedef const float CFloatingTypeGPU;
//typedef double FloatingTypeGPU;

#define USE_GPU_COLORING_DATA 
//#define EVAL_TET_ACCESS_COUNT 

namespace GAIA {
	// a struct that feeds to GPU,
	// contains pointers to TetMesh's GPU data buffers

	struct TetMeshTopologyGPU final {
		//   tet's vertex buffer: 4 x nTets
		int32_t* tetVIds;

		int32_t* vertexNeighborTets;
		int32_t* vertexNeighborTets_vertexOrder;
		int32_t* vertexNeighborTets_infos;

		int32_t nVerts;
		int32_t nTets;
	};

	struct TetMeshFEMGPU 
	{
		typedef TetMeshFEMGPU* Ptr;

		// topology
		TetMeshTopologyGPU* pTopology;

		//* non-variant data
		// data arranged per vertex
		FloatingTypeGPU* vertexMass;

		// data arranged per tet
		FloatingTypeGPU* DmInvs;            // tet's mats, 9 x nTets
		FloatingTypeGPU* tetRestVolume;
		FloatingTypeGPU* tetInvRestVolume;
		//* non-variant data

		// flatten vertex data: 3*nVerts
		FloatingTypeGPU* vertPos;
		FloatingTypeGPU* vertPrevPos;
		FloatingTypeGPU* velocity;

		int8_t* vertexFixedMask;

		// // tet's inversion sign buffer: nTets
		// int8_t* tetInvertedSign;
		// int8_t* verticesInvertedSign;

#ifdef EVAL_TET_ACCESS_COUNT
		int32_t* tetAccessCount;
		int32_t* vertAccessCount;
#endif // !EVAL_TET_ACCESS_COUNT

		bool activeForCollision = true;
		bool activeForMaterialSolve = true;

	};

	struct PhysicsDataGPU
	{
		// non-variant data
		FloatingTypeGPU dt;

	};

	GPU_CPU_INLINE_FUNC int32_t getVertexNeighborTetsStart(TetMeshTopologyGPU* pToplogy, int32_t vertexId) {
		return pToplogy->vertexNeighborTets[2 * vertexId];
	}

	GPU_CPU_INLINE_FUNC int32_t getVertexNeighborTetsNumber(TetMeshTopologyGPU* pToplogy, int32_t vertexId) {
		return pToplogy->vertexNeighborTets[2 * vertexId + 1];
	}

}
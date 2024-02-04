#pragma once

#define VERTEX_BUFFER_STRIDE 3

#include "CuMatrix/Buffers/ManagedBuffer.h"
typedef float FloatingTypeGPU;
//typedef double FloatingTypeGPU;

#define USE_GPU_COLORING_DATA 
//#define EVAL_TET_ACCESS_COUNT 

namespace GAIA {
	// a struct that feeds to GPU,
	// contains pointers to TetMesh's GPU data buffers
	struct PBDTetMeshFEMGPU 
	{
		typedef std::shared_ptr<PBDTetMeshFEMGPU> SharedPtr;
		typedef PBDTetMeshFEMGPU* Ptr;

		// non-variant data
		FloatingTypeGPU* tetRestVolume;
		FloatingTypeGPU* tetInvRestVolume;
		FloatingTypeGPU* vertexInvMass;

		//   tet's mats, 9 x nTets
		FloatingTypeGPU* DmInvs;
		//   tet's vertex buffer: 4 x nTets
		const int32_t* tetVIds;
		// tet's inversion sign buffer: nTets
		int8_t* tetInvertedSign;
		int8_t* verticesInvertedSign;

		// flatten vertex data: 3 nVerts
		FloatingTypeGPU* vertPos;
		FloatingTypeGPU* vertPrevPos;
		FloatingTypeGPU* velocity;

		// non-variant data
		FloatingTypeGPU dt;
		int32_t nVerts;
		int32_t nTets;


#ifdef EVAL_TET_ACCESS_COUNT
		int32_t* tetAccessCount;
		int32_t* vertAccessCount;
#endif // !USE_GPU_COLORING_DATA

#ifdef USE_GPU_COLORING_DATA
		// tet coloring
		// those information can be feed by the CPU when initializing a parallelization group
		int32_t numTetsColoredCatergories;
		int32_t* tetsColoringEachCategorySize;
		int32_t** tetsColoringCategoriesPointers;
#endif // !USE_GPU_COLORING_DATA


	};

}
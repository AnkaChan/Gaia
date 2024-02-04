#pragma once
#include <stdint.h>
#include "PBDTetMeshFEMGPU.h"

namespace GAIA {
	template <class Func>
	__global__ void parallel_for_tets(int32_t* tetsParallelizationGroup, int32_t numTets, PBDTetMeshFEMGPU* pTetMeshGPU, Func func) {
		for (int iTet = blockIdx.x * blockDim.x + threadIdx.x;
			iTet < numTets;
			iTet += blockDim.x * gridDim.x)
		{
			int32_t tetId = tetsParallelizationGroup[iTet];
			func(tetId, pTetMeshGPU);
		}
	}

	template <class Func>
	__global__ void parallel_for_tets(PBDTetMeshFEMGPU* pTetMeshGPU, Func func) {
		for (int iTet = blockIdx.x * blockDim.x + threadIdx.x;
			iTet < pTetMeshGPU->nTets;
			iTet += blockDim.x * gridDim.x)
		{
			func(iTet, pTetMeshGPU);
		}
	}


	template <class Func>
	__global__ void parallel_for_verts(int32_t numVerts, PBDTetMeshFEMGPU* pTetMeshGPU, Func func) {
		for (int iVert = blockIdx.x * blockDim.x + threadIdx.x;
			iVert < numVerts;
			iVert += blockDim.x * gridDim.x)
		{
			func(iVert, pTetMeshGPU);
		}
	}
}
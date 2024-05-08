#include "VBD_NeoHookeanGPU.h"
#include "VBDPhysicsCompute.h"
#include "VBD_GeneralCompute.h"
#include "VBDPhysicsTest.h"

#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrixVis.h"
#include "../../3rdParty/CuMatrix/CuMatrix/Geometry/Geometry.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"  
#define GRAVITY_DIM 1
//#define PRINT_DBG_INFO
//#define USE_IPC_FRICTION

#define VEC12_SET(vecOut, vecIn)\
vecOut[0]=vecIn[0]; vecOut[1]=vecIn[1]; vecOut[2]=vecIn[2];\
vecOut[3]=vecIn[3]; vecOut[4]=vecIn[4]; vecOut[5]=vecIn[5];\
vecOut[6]=vecIn[6]; vecOut[7]=vecIn[7]; vecOut[8]=vecIn[8];\
vecOut[9]=vecIn[9]; vecOut[10]=vecIn[10]; vecOut[11]=vecIn[11];

#define VEC12_SET_CONSTANT(vecOut, c)\
vecOut[0]=c; vecOut[1]=c; vecOut[2]=c;\
vecOut[3]=c; vecOut[4]=c; vecOut[5]=c;\
vecOut[6]=c; vecOut[7]=c; vecOut[8]=c;\
vecOut[9]=c; vecOut[10]=c; vecOut[11]=c;

#define VEC12_ADD(a, b, out)\
out[0]=a[0] + b[0]; out[1]=a[1] +   b[1];  out[2]=a[2] + b[2];\
out[3]=a[3] + b[3]; out[4]=a[4] +   b[4];  out[5]=a[5] + b[5];\
out[6]=a[6] + b[6]; out[7]=a[7] +   b[7];  out[8]=a[8] + b[8];\
out[9]=a[9] + b[9]; out[10]=a[10] + b[10]; out[11]=a[11] + b[11];


using namespace GAIA;

// compute vertex force and hessian based on gradient and Hessian of F, and accumulate them to force and h
template<template<class> class MatType>
GPU_CPU_INLINE_FUNC  void assembleVertexVForceAndHessian(CFloatingTypeGPU* dE_dF, MatType<FloatingTypeGPU>* d2E_dF_dF, CFloatingTypeGPU m1, CFloatingTypeGPU m2, CFloatingTypeGPU m3,
	FloatingTypeGPU* force, FloatingTypeGPU* h)
{
	CFloatingTypeGPU A1 = dE_dF[0];
	CFloatingTypeGPU A2 = dE_dF[1];
	CFloatingTypeGPU A3 = dE_dF[2];
	CFloatingTypeGPU A4 = dE_dF[3];
	CFloatingTypeGPU A5 = dE_dF[4];
	CFloatingTypeGPU A6 = dE_dF[5];
	CFloatingTypeGPU A7 = dE_dF[6];
	CFloatingTypeGPU A8 = dE_dF[7];
	CFloatingTypeGPU A9 = dE_dF[8];

	// force is the negative of gradient
	force[0] -= A1 * m1 + A4 * m2 + A7 * m3;
	force[1] -= A2 * m1 + A5 * m2 + A8 * m3;
	force[2] -= A3 * m1 + A6 * m2 + A9 * m3;

	/*const CuMatrix::Mat9x9Abstract<FloatingTypeGPU>& H = *d2E_dF_dF;
	h[0] += m1 * (H(0, 0) * m1 + H(3, 0) * m2 + H(6, 0) * m3) + m2 * (H(0, 3) * m1 + H(3, 3) * m2 + H(6, 3) * m3) + m3 * (H(0, 6) * m1 + H(3, 6) * m2 + H(6, 6) * m3);
	h[1] += m1 * (H(1, 0) * m1 + H(4, 0) * m2 + H(7, 0) * m3) + m2 * (H(1, 3) * m1 + H(4, 3) * m2 + H(7, 3) * m3) + m3 * (H(1, 6) * m1 + H(4, 6) * m2 + H(7, 6) * m3);
	h[2] += m1 * (H(2, 0) * m1 + H(5, 0) * m2 + H(8, 0) * m3) + m2 * (H(2, 3) * m1 + H(5, 3) * m2 + H(8, 3) * m3) + m3 * (H(2, 6) * m1 + H(5, 6) * m2 + H(8, 6) * m3);
	h[3] += m1 * (H(0, 1) * m1 + H(3, 1) * m2 + H(6, 1) * m3) + m2 * (H(0, 4) * m1 + H(3, 4) * m2 + H(6, 4) * m3) + m3 * (H(0, 7) * m1 + H(3, 7) * m2 + H(6, 7) * m3);
	h[4] += m1 * (H(1, 1) * m1 + H(4, 1) * m2 + H(7, 1) * m3) + m2 * (H(1, 4) * m1 + H(4, 4) * m2 + H(7, 4) * m3) + m3 * (H(1, 7) * m1 + H(4, 7) * m2 + H(7, 7) * m3);
	h[5] += m1 * (H(2, 1) * m1 + H(5, 1) * m2 + H(8, 1) * m3) + m2 * (H(2, 4) * m1 + H(5, 4) * m2 + H(8, 4) * m3) + m3 * (H(2, 7) * m1 + H(5, 7) * m2 + H(8, 7) * m3);
	h[6] += m1 * (H(0, 2) * m1 + H(3, 2) * m2 + H(6, 2) * m3) + m2 * (H(0, 5) * m1 + H(3, 5) * m2 + H(6, 5) * m3) + m3 * (H(0, 8) * m1 + H(3, 8) * m2 + H(6, 8) * m3);
	h[7] += m1 * (H(1, 2) * m1 + H(4, 2) * m2 + H(7, 2) * m3) + m2 * (H(1, 5) * m1 + H(4, 5) * m2 + H(7, 5) * m3) + m3 * (H(1, 8) * m1 + H(4, 8) * m2 + H(7, 8) * m3);
	h[8] += m1 * (H(2, 2) * m1 + H(5, 2) * m2 + H(8, 2) * m3) + m2 * (H(2, 5) * m1 + H(5, 5) * m2 + H(8, 5) * m3) + m3 * (H(2, 8) * m1 + H(5, 8) * m2 + H(8, 8) * m3);*/

	FloatingTypeGPU HL_Row1[9];
	FloatingTypeGPU HL_Row2[9];
	FloatingTypeGPU HL_Row3[9];

	FloatingTypeGPU* HL[3] = { HL_Row1, HL_Row2, HL_Row3 };

	for (int32_t iCol = 0; iCol < 9; iCol++)
	{
		HL_Row1[iCol] = (*d2E_dF_dF)(0, iCol) * m1;
		HL_Row1[iCol] += (*d2E_dF_dF)(3, iCol) * m2;
		HL_Row1[iCol] += (*d2E_dF_dF)(6, iCol) * m3;

		HL_Row2[iCol] = (*d2E_dF_dF)(1, iCol) * m1;
		HL_Row2[iCol] += (*d2E_dF_dF)(4, iCol) * m2;
		HL_Row2[iCol] += (*d2E_dF_dF)(7, iCol) * m3;


		HL_Row3[iCol] = (*d2E_dF_dF)(2, iCol) * m1;
		HL_Row3[iCol] += (*d2E_dF_dF)(5, iCol) * m2;
		HL_Row3[iCol] += (*d2E_dF_dF)(8, iCol) * m3;
	}

	for (int32_t iRow = 0; iRow < 3; iRow++) {
		h[iRow] += HL[iRow][0] * m1 + HL[iRow][3] * m2 + HL[iRow][6] * m3;
		h[iRow + 3] += HL[iRow][1] * m1 + HL[iRow][4] * m2 + HL[iRow][7] * m3;
		h[iRow + 6] += HL[iRow][2] * m1 + HL[iRow][5] * m2 + HL[iRow][8] * m3;
	}

}

void GAIA::VBDSolveParallelGroup_collisionSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t numActiveCollisions,
	int32_t* activeCollisionsEachParallelGroup, int32_t numThreads, cudaStream_t cudaStream)
{
	VBDSolveParallelGroup_collisionSweep_kernel KERNEL_ARGS4((numActiveCollisions + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, numActiveCollisions, activeCollisionsEachParallelGroup);
}

void GAIA::VBDSolveParallelGroup_tetSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead,
	int32_t tetParallelGroupSize, int32_t numThreads, cudaStream_t cudaStream)
{
	// VBDSolveParallelGroup_tetSweep_kernel KERNEL_ARGS4((tetParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
	// 	(pPhysicsData, tetParallelGroupsHead, tetParallelGroupSize);

	//std::cout << "tetSweep_kernel: " << (tetParallelGroupSize + NUM_THREADS_TET_SWEEP - 1) / NUM_THREADS_TET_SWEEP << " blocks" << NUM_THREADS_TET_SWEEP << " threads\n";

	VBDSolveParallelGroup_tetSweep_kernel_V2 KERNEL_ARGS4((tetParallelGroupSize + NUM_THREADS_TET_SWEEP - 1) / NUM_THREADS_TET_SWEEP, NUM_THREADS_TET_SWEEP, 0, cudaStream)
		(pPhysicsData, tetParallelGroupsHead, tetParallelGroupSize);
}

void GAIA::VBDSolveParallelGroup_vertexSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead,
	int32_t veretxParallelGroupSize, CFloatingTypeGPU acceleratorOmega, int32_t numThreads, cudaStream_t cudaStream)
{
	//VBDSolveParallelGroup_vertexSweep_kernel KERNEL_ARGS4((veretxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
	//	(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize, acceleratorOmega);

	// V2, one block for each vertex
	//std::cout << "vertexSweep_kernel: " << (veretxParallelGroupSize + NUM_THREADS_TET_SWEEP - 1) / NUM_THREADS_TET_SWEEP << " blocks" << NUM_THREADS_TET_SWEEP << " threads\n";
	//VBDSolveParallelGroup_vertexSweep_kernel_V2 KERNEL_ARGS4(veretxParallelGroupSize, NUM_THREADS_VERTEX_SWEEP, 0, cudaStream)
	//	(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize, acceleratorOmega);

	VBDSolveParallelGroup_allInOne_kernel_V2 KERNEL_ARGS4(veretxParallelGroupSize, NUM_THREADS_VERTEX_SWEEP, 0, cudaStream)
		(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize, acceleratorOmega);
}

void GAIA::VBDSolveParallelGroup_vertexSweepAcceleratedGS(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize, FloatingTypeGPU acceleratorOmega, int32_t numThreads, cudaStream_t cudaStream)
{
	VBDSolveParallelGroup_vertexSweepAcceleratedGS_kernel KERNEL_ARGS4((veretxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize, acceleratorOmega);
}

void GAIA::VBDSolveParallelGroup_allInOneSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, FloatingTypeGPU acceleratorOmega,
	int32_t numThreads, cudaStream_t cudaStream)
{
	VBDSolveParallelGroup_allInOne_kernel KERNEL_ARGS4((parallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, parallelGroupsHead, parallelGroupSize, acceleratorOmega);

	// VBDSolveParallelGroup_kernel_serial KERNEL_ARGS2(1, 1)	(pPhysicsData, parallelGroupsHead, parallelGroupSize);
}

void GAIA::GDSolveParallelGroup_allInOneSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize, 
	int32_t numThreads, cudaStream_t cudaStream)
{
	GDSolveParallelGroup_vertexSweep_kernel KERNEL_ARGS4((veretxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize);
}

void GAIA::GDSolveParallelGroup_BlockJacobi_allInOneSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize, int32_t numThreads, cudaStream_t cudaStream)
{
	GDSolveParallelGroup_BlockJacobi_allInOneSweep_kernel KERNEL_ARGS4(veretxParallelGroupSize, NUM_THREADS_VERTEX_SWEEP, 0, cudaStream)
		(pPhysicsData, vertexParallelGroupsHead, veretxParallelGroupSize);
}

void GAIA::GDSolveParallelGroup_updatePositionSweep(VBDPhysicsDataGPU* pPhysicsData, CFloatingTypeGPU stepSize, CFloatingTypeGPU acceleratorOmega, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize, int32_t numThreads, cudaStream_t cudaStream)
{
	GDSolveParallelGroup_updatePositionSweep_kernel KERNEL_ARGS4((veretxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, stepSize, acceleratorOmega, vertexParallelGroupsHead, veretxParallelGroupSize);
}


void GAIA::VBDUpdateVelocity(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices,
	int32_t numThreads, cudaStream_t cudaStream)
{
	VBDUpdateRelativeVelocity_kernel KERNEL_ARGS4((numVertices + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, vertexAllParallelGroups, numVertices);

	VBDUpdateVelocity_kernel KERNEL_ARGS4((numVertices + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, vertexAllParallelGroups, numVertices);
}

void GAIA::evaluateElasticEnergyGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize, FloatingTypeGPU* tetElasticEnergy, int32_t numThreads, cudaStream_t cudaStream)
{
	evaluateElasticEnergy_kernel KERNEL_ARGS4((tetParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, cudaStream)
		(pPhysicsData, tetParallelGroupsHead, tetParallelGroupSize, tetElasticEnergy);
}

__global__ void GAIA::VBDUpdateVelocity_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < numVertices)
	{
		int32_t meshId = vertexAllParallelGroups[iVert * 2];
		int32_t vertexId = vertexAllParallelGroups[iVert * 2 + 1];
		VBDBaseTetMeshGPU* pTetMeshGPU = (VBDBaseTetMeshGPU*)pPhysicsData->tetMeshes[meshId];
		FloatingTypeGPU* vel = pTetMeshGPU->velocity + VERTEX_BUFFER_STRIDE * vertexId;

		if (pTetMeshGPU->activeForMaterialSolve)
		{
			if (pTetMeshGPU->vertexFixedMask[vertexId])
			{
				CuMatrix::vec3Set(vel, 0.f);
			}
			else {
				CFloatingTypeGPU dt = pPhysicsData->dt;
				CFloatingTypeGPU vertInvMass = 1.f / pTetMeshGPU->vertexMass[vertexId];
				updateVertexVelocity(pTetMeshGPU, vertexId, dt);
#ifndef USE_IPC_BOUNDARY_FRICTION
				applyBoundaryFriction(pPhysicsData, pTetMeshGPU, vertexId, dt, vertInvMass);
#endif // !USE_IPC_BOUNDARY_FRICTION
#ifndef USE_IPC_FRICTION
				applyCollisionFriction(pPhysicsData, pTetMeshGPU, vertexId, dt, vertInvMass);
#endif // !USE_IPC_FRICTION
				dampVelocity(vel, pTetMeshGPU->maxVelocityMagnitude, pTetMeshGPU->exponentialVelDamping, pTetMeshGPU->constantVelDamping);
			}
		}

	}
	return;
}

__global__ void GAIA::VBDApplyAccelerator_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices, FloatingTypeGPU acceleratorOmega)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < numVertices)
	{
		int32_t meshId = vertexAllParallelGroups[iVert * 2];
		int32_t vertexId = vertexAllParallelGroups[iVert * 2 + 1];


	}
}

__global__ void GAIA::VBDUpdateRelativeVelocity_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < numVertices)
	{
		int32_t meshId = vertexAllParallelGroups[iVert * 2];
		int32_t vertexId = vertexAllParallelGroups[iVert * 2 + 1];
		computeRelativeVelocityGPU(pPhysicsData, meshId, vertexId);
	}
	return;
}

__global__ void GAIA::VBDSolveParallelGroup_allInOne_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < parallelGroupSize)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];

		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

		if (!pTetMeshGPU->vertexFixedMask[vertexId]
			&& pTetMeshGPU->activeForMaterialSolve
			)
		{
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);
			// compute y_k+1
			VBDStep_allInOne(pPhysicsData, pTetMeshGPU, vertexId);

			if (pPhysicsData->useAccelerator
				&& acceleratorOmega != 1.0f)
			{
				// compute y_k+1

				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;

				// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

				FloatingTypeGPU* y_new = pTetMeshGPU->positionsNew + VERTEX_BUFFER_STRIDE * vertexId;

				FloatingTypeGPU yDiff[3];
				CuMatrix::vec3Minus(v, y_prev, yDiff);
				CuMatrix::vec3Set(y_new, y_prev);
				CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, y_new);

				// write y_k to y_k-1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}

	}
	return;
}

__global__ void GAIA::VBDSolveParallelGroup_collisionSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t numActiveCollisions,
	int32_t* activeCollisionsEachParallelGroup)
{
	int32_t iActiveCollision = blockIdx.x * blockDim.x + threadIdx.x;
	if (iActiveCollision < numActiveCollisions)
	{

		int32_t meshId = activeCollisionsEachParallelGroup[iActiveCollision * 2];
		int32_t vertexId = activeCollisionsEachParallelGroup[iActiveCollision * 2 + 1];

		if (pPhysicsData->tetMeshes[meshId]->activeForCollision)
		{
			updateCollisionInfoGPU(pPhysicsData, meshId, vertexId);
		}
	}

	return;
}

__global__ void GAIA::VBDSolveParallelGroup_tetSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize)
{
	int32_t iTet = blockIdx.x * blockDim.x + threadIdx.x;
	if (iTet < tetParallelGroupSize)
	{
		int32_t meshId = tetParallelGroupsHead[iTet * 4];
		int32_t tetId = tetParallelGroupsHead[iTet * 4 + 1];
		int32_t vertexOrderInTet = tetParallelGroupsHead[iTet * 4 + 2];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		int32_t* tetVIds = pTetMeshGPU->pTopology->tetVIds + 4 * tetId;

		int32_t vId = tetVIds[vertexOrderInTet];
		CFloatingTypeGPU* v = pTetMeshGPU->getVert(vId);

		if (pTetMeshGPU->activeForMaterialSolve)
		{

			FloatingTypeGPU tetFAndH[12] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };

			evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[tetId], tetId, vertexOrderInTet,
				tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, tetFAndH, tetFAndH + 3, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysicsData->dt);

			FloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * tetId;

			for (int32_t i = 0; i < 12; i++)
			{
				tetFAndHStorage[i] = tetFAndH[i];
			}
		}
	}

	return;
}

__global__ void GAIA::VBDSolveParallelGroup_tetSweep_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize)
{
	__shared__  __builtin_align__(16) FloatingTypeGPU tetFAndHAll[12 * NUM_THREADS_TET_SWEEP];
	//__shared__  CuMatrix::Mat9x9Static<FloatingTypeGPU> d2E_dF_dF_all[NUM_THREADS_TET_SWEEP];

	
	int32_t iTet = blockIdx.x * blockDim.x + threadIdx.x;
	int32_t tid = threadIdx.x;

	if (iTet < tetParallelGroupSize)
	{
		int32_t meshId = tetParallelGroupsHead[iTet * 4];
		int32_t tetId = tetParallelGroupsHead[iTet * 4 + 1];
		int32_t vertexOrderInTet = tetParallelGroupsHead[iTet * 4 + 2];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		int32_t* tetVIds = pTetMeshGPU->pTopology->tetVIds + 4 * tetId;

		int32_t vId = tetVIds[vertexOrderInTet];
		CFloatingTypeGPU* v = pTetMeshGPU->getVert(vId);

		if (pTetMeshGPU->activeForMaterialSolve)
		{
			CFloatingTypeGPU* vs = pTetMeshGPU->vertPos;
			CFloatingTypeGPU* vsPrev = pTetMeshGPU->vertPrevPos;
			CFloatingTypeGPU miu = pTetMeshGPU->miu;
			CFloatingTypeGPU lmbd = pTetMeshGPU->lmbd;

			FloatingTypeGPU* tetFAndH = tetFAndHAll + tid * 12;
			//FloatingTypeGPU tetFAndH[12];
			VEC12_SET_CONSTANT(tetFAndH, 0.f);
			FloatingTypeGPU* force = tetFAndH;
			FloatingTypeGPU* h = tetFAndH + 3;

			FloatingTypeGPU dampingVolume = pTetMeshGPU->dampingVolume;
			FloatingTypeGPU dampingShear = pTetMeshGPU->dampingShear;

			CFloatingTypeGPU* DmInvs = pTetMeshGPU->DmInvs;

			CFloatingTypeGPU A = pTetMeshGPU->tetRestVolume[tetId];
			CFloatingTypeGPU dt = pPhysicsData->dt;

			CFloatingTypeGPU* v0 = vs + VERTEX_BUFFER_STRIDE * tetVIds[0];
			CFloatingTypeGPU* v1 = vs + VERTEX_BUFFER_STRIDE * tetVIds[1];
			CFloatingTypeGPU* v2 = vs + VERTEX_BUFFER_STRIDE * tetVIds[2];
			CFloatingTypeGPU* v3 = vs + VERTEX_BUFFER_STRIDE * tetVIds[3];

			FloatingTypeGPU Ds[9];
			FloatingTypeGPU F[9];

			CuMatrix::vec3Minus(v1, v0, Ds);
			CuMatrix::vec3Minus(v2, v0, Ds + 3);
			CuMatrix::vec3Minus(v3, v0, Ds + 6);
			CFloatingTypeGPU* DmInv = DmInvs + tetId * 9;
			//printf("Ds:\n");
			//CuMatrix::printMat3(Ds);

			CuMatrix::mat3MatProduct(Ds, DmInv, F);
			//printf("F:\n");
			//CuMatrix::printMat3(F);

			CFloatingTypeGPU a = 1 + miu / lmbd;

			//CFloatingTypeGPU Phi_D = 0.5f * (CuMatrix::mat3FNormSquare(F) - 3.f);
			CFloatingTypeGPU detF = CuMatrix::mat3Determinant(F);
			//CFloatingTypeGPU Phi_H = 0.5f * SQR(detF - a);

			// CFloatingTypeGPU E = A * (miu * Phi_D + lmbd * Phi_H);

			// printf("Phi_D: %f | Phi_H: %f | E: %f\n", Phi_D, detF, Phi_H);

			FloatingTypeGPU* dPhi_D_dF = F;

			// printf("dPhi_D_dF:\n");
			// CuMatrix::printFloatVec(dPhi_D_dF, 9);

			CFloatingTypeGPU F1_1 = F[0];
			CFloatingTypeGPU F2_1 = F[1];
			CFloatingTypeGPU F3_1 = F[2];
			CFloatingTypeGPU F1_2 = F[3 * 1];
			CFloatingTypeGPU F2_2 = F[1 + 3 * 1];
			CFloatingTypeGPU F3_2 = F[2 + 3 * 1];
			CFloatingTypeGPU F1_3 = F[3 * 2];
			CFloatingTypeGPU F2_3 = F[1 + 3 * 2];
			CFloatingTypeGPU F3_3 = F[2 + 3 * 2];

			FloatingTypeGPU ddetF_dF[9] =
			{ F2_2 * F3_3 - F2_3 * F3_2,
			  F1_3 * F3_2 - F1_2 * F3_3,
			  F1_2 * F2_3 - F1_3 * F2_2,
			  F2_3 * F3_1 - F2_1 * F3_3,
			  F1_1 * F3_3 - F1_3 * F3_1,
			  F1_3 * F2_1 - F1_1 * F2_3,
			  F2_1 * F3_2 - F2_2 * F3_1,
			  F1_2 * F3_1 - F1_1 * F3_2,
			  F1_1 * F2_2 - F1_2 * F2_1 };

			// printf("ddetF_dF:\n");
			// CuMatrix::printFloatVec(ddetF_dF, 9);

			//CuMatrix::Mat9x9Static<FloatingTypeGPU>& d2E_dF_dF = d2E_dF_dF_all[tid];
			CuMatrix::Mat9x9Static<FloatingTypeGPU> d2E_dF_dF;

			CuMatrix::vec9OuterProduct(ddetF_dF, ddetF_dF, d2E_dF_dF);

			CFloatingTypeGPU k = detF - a;
			d2E_dF_dF(0, 4) += k * F3_3;
			d2E_dF_dF(4, 0) += k * F3_3;
			d2E_dF_dF(0, 5) += k * -F2_3;
			d2E_dF_dF(5, 0) += k * -F2_3;
			d2E_dF_dF(0, 7) += k * -F3_2;
			d2E_dF_dF(7, 0) += k * -F3_2;
			d2E_dF_dF(0, 8) += k * F2_2;
			d2E_dF_dF(8, 0) += k * F2_2;

			d2E_dF_dF(1, 3) += k * -F3_3;
			d2E_dF_dF(3, 1) += k * -F3_3;
			d2E_dF_dF(1, 5) += k * F1_3;
			d2E_dF_dF(5, 1) += k * F1_3;
			d2E_dF_dF(1, 6) += k * F3_2;
			d2E_dF_dF(6, 1) += k * F3_2;
			d2E_dF_dF(1, 8) += k * -F1_2;
			d2E_dF_dF(8, 1) += k * -F1_2;

			d2E_dF_dF(2, 3) += k * F2_3;
			d2E_dF_dF(3, 2) += k * F2_3;
			d2E_dF_dF(2, 4) += k * -F1_3;
			d2E_dF_dF(4, 2) += k * -F1_3;
			d2E_dF_dF(2, 6) += k * -F2_2;
			d2E_dF_dF(6, 2) += k * -F2_2;
			d2E_dF_dF(2, 7) += k * F1_2;
			d2E_dF_dF(7, 2) += k * F1_2;

			d2E_dF_dF(3, 7) += k * F3_1;
			d2E_dF_dF(7, 3) += k * F3_1;
			d2E_dF_dF(3, 8) += k * -F2_1;
			d2E_dF_dF(8, 3) += k * -F2_1;

			d2E_dF_dF(4, 6) += k * -F3_1;
			d2E_dF_dF(6, 4) += k * -F3_1;
			d2E_dF_dF(4, 8) += k * F1_1;
			d2E_dF_dF(8, 4) += k * F1_1;

			d2E_dF_dF(5, 6) += k * F2_1;
			d2E_dF_dF(6, 5) += k * F2_1;
			d2E_dF_dF(5, 7) += k * -F1_1;
			d2E_dF_dF(7, 5) += k * -F1_1;

			d2E_dF_dF.multiplyBy(lmbd);

			d2E_dF_dF(0, 0) += miu;
			d2E_dF_dF(1, 1) += miu;
			d2E_dF_dF(2, 2) += miu;
			d2E_dF_dF(3, 3) += miu;
			d2E_dF_dF(4, 4) += miu;
			d2E_dF_dF(5, 5) += miu;
			d2E_dF_dF(6, 6) += miu;
			d2E_dF_dF(7, 7) += miu;
			d2E_dF_dF(8, 8) += miu;

			d2E_dF_dF.multiplyBy(A);

			// printf("d2E_dF_dF:\n");
			// CuMatrix::printMat(d2E_dF_dF.data, 9, 9);

			//std::cout << "d2E_dF_dF:\n" << d2E_dF_dF << std::endl;

			CuMatrix::vec9Mul(dPhi_D_dF, miu, dPhi_D_dF);
			CuMatrix::vec9Mul(ddetF_dF, lmbd * k, ddetF_dF);

			FloatingTypeGPU dE_dF[9];

			CuMatrix::vec9Add(dPhi_D_dF, ddetF_dF, dE_dF);

			CuMatrix::vec9Mul(dE_dF, A, dE_dF);

			/*printf("dE_dF of tet %d:\n", tetId);
			CuMatrix::printFloatVec(dE_dF, 9);*/

			//std::cout << "dE_dF:\n" << dE_dF << std::endl;

			CFloatingTypeGPU DmInv1_1 = DmInv[0];
			CFloatingTypeGPU DmInv2_1 = DmInv[1];
			CFloatingTypeGPU DmInv3_1 = DmInv[2];
			CFloatingTypeGPU DmInv1_2 = DmInv[3 * 1];
			CFloatingTypeGPU DmInv2_2 = DmInv[1 + 3 * 1];
			CFloatingTypeGPU DmInv3_2 = DmInv[2 + 3 * 1];
			CFloatingTypeGPU DmInv1_3 = DmInv[3 * 2];
			CFloatingTypeGPU DmInv2_3 = DmInv[1 + 3 * 2];
			CFloatingTypeGPU DmInv3_3 = DmInv[2 + 3 * 2];

			CFloatingTypeGPU ms[4][3] = {
				{-DmInv1_1 - DmInv2_1 - DmInv3_1, -DmInv1_2 - DmInv2_2 - DmInv3_2,  -DmInv1_3 - DmInv2_3 - DmInv3_3},
				{DmInv1_1, DmInv1_2, DmInv1_3},
				{DmInv2_1, DmInv2_2, DmInv2_3},
				{DmInv3_1, DmInv3_2, DmInv3_3}
			};
			CFloatingTypeGPU m1 = ms[vertexOrderInTet][0], m2 = ms[vertexOrderInTet][1], m3 = ms[vertexOrderInTet][2];
			assembleVertexVForceAndHessian(dE_dF, &d2E_dF_dF, m1, m2, m3, force, h);

			if (dampingVolume > 0.f || dampingShear > 0.f)
			{
				FloatingTypeGPU dampingH[9];
				CuMatrix::vec9Mul(h, dampingVolume, dampingH);

				FloatingTypeGPU tmp = (m1 * m1 + m2 * m2 + m3 * m3) * miu * A;

				h[0] += tmp;
				h[4] += tmp;
				h[8] += tmp;
				tmp *= dampingShear;
				dampingH[0] += tmp;
				dampingH[4] += tmp;
				dampingH[8] += tmp;
				CuMatrix::vec9Mul(dampingH, 1 / dt, dampingH);

				CFloatingTypeGPU* v = vs + VERTEX_BUFFER_STRIDE * tetVIds[vertexOrderInTet];
				CFloatingTypeGPU* vPrev = vsPrev + VERTEX_BUFFER_STRIDE * tetVIds[vertexOrderInTet];
				FloatingTypeGPU displacement[3];
				CuMatrix::vec3Minus(vPrev, v, displacement);

				FloatingTypeGPU dampingForce[3];
				CuMatrix::mat3VecProduct(dampingH, displacement, dampingForce);
				CuMatrix::vec3Add(force, dampingForce, force);

				CuMatrix::vec9Add(h, dampingH, h);
			}

			FloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * tetId;

			for (int32_t i = 0; i < 12; i++)
			{
				tetFAndHStorage[i] = tetFAndH[i];
			}
		}
	}

	return;
}

__global__ void GAIA::GDSolveParallelGroup_tetSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize)
{
	int32_t iTet = blockIdx.x * blockDim.x + threadIdx.x;
	if (iTet < tetParallelGroupSize)
	{
		int32_t meshId = tetParallelGroupsHead[iTet * 4];
		int32_t tetId = tetParallelGroupsHead[iTet * 4 + 1];
		int32_t vertexOrderInTet = tetParallelGroupsHead[iTet * 4 + 2];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		int32_t* tetVIds = pTetMeshGPU->pTopology->tetVIds + 4 * tetId;

		int32_t vId = tetVIds[vertexOrderInTet];
		CFloatingTypeGPU* v = pTetMeshGPU->getVert(vId);

		if (pTetMeshGPU->activeForMaterialSolve)
		{

			FloatingTypeGPU tetFAndH[6] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};

			evaluateNeoHookeanMaterialForceAndDiagonalHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[tetId], tetId, vertexOrderInTet,
				tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->DmInvs, tetFAndH, tetFAndH + 3);

			FloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * tetId;

			for (int32_t i = 0; i < 6; i++)
			{
				tetFAndHStorage[i] = tetFAndH[i];
			}
		}
	}

	return;
}


__global__ void GAIA::VBDSolveParallelGroup_vertexSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < parallelGroupSize)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

		if (!pTetMeshGPU->vertexFixedMask[vertexId]
			&& pTetMeshGPU->activeForMaterialSolve
			//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
			)
		{
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);
			// compute y_k+1
			VBDStep_vertexSweep(pPhysicsData, pTetMeshGPU, vertexId, acceleratorOmega);

			if (pPhysicsData->useAccelerator
				&& acceleratorOmega != 1.0f)
			{
				// compute y_k+1
				FloatingTypeGPU* y_new = pTetMeshGPU->positionsNew + VERTEX_BUFFER_STRIDE * vertexId;
				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;
				if (pTetMeshGPU->activeCollisionMask[vertexId])
				{
					CuMatrix::vec3Set(y_new, v);
				}
				else
				{
					// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

					FloatingTypeGPU yDiff[3];
					CuMatrix::vec3Minus(v, y_prev, yDiff);
					CuMatrix::vec3Set(y_new, y_prev);
					CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, y_new);

				}
				// write y_k to y_k-1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}

	}
	return;
}

__inline__ __device__ void unwarppedReduction(volatile FloatingTypeGPU* forceAndHessianAll, int32_t tid)
{
	// for NUM_THREADS_VERTEX_SWEEP == 16 only

	int32_t index = tid * 12;
	volatile FloatingTypeGPU* forceAndHessian = forceAndHessianAll + index;

	VEC12_ADD(forceAndHessian, (forceAndHessian + 8 * 12), forceAndHessian);
	VEC12_ADD(forceAndHessian, (forceAndHessian + 4 * 12), forceAndHessian);
	VEC12_ADD(forceAndHessian, (forceAndHessian + 2 * 12), forceAndHessian);
	VEC12_ADD(forceAndHessian, (forceAndHessian + 1 * 12), forceAndHessian);

	//forceAndHessianAll[tid] += s_data[tid + 32];
	//forceAndHessianAll[tid] += s_data[tid + 16];
	//forceAndHessianAll[tid] += s_data[tid + 8];
	//forceAndHessianAll[tid] += s_data[tid + 4];
	//forceAndHessianAll[tid] += s_data[tid + 2];
	//forceAndHessianAll[tid] += s_data[tid + 1];
}

__global__ void GAIA::VBDSolveParallelGroup_vertexSweep_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega)
{
	__shared__ __builtin_align__(16) FloatingTypeGPU forceAndHessianAll[NUM_THREADS_VERTEX_SWEEP * 12];

	int32_t iV = blockIdx.x;
	int32_t meshId = parallelGroupsHead[iV * 2];
	int32_t vertexId = parallelGroupsHead[iV * 2 + 1];
	VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
	FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

	if (!pTetMeshGPU->vertexFixedMask[vertexId]
		&& pTetMeshGPU->activeForMaterialSolve
		//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
		)
	{
		int32_t tid = threadIdx.x;

		// assign material force & hessian to shared memory

		FloatingTypeGPU* forceAndHessian = forceAndHessianAll + 12 * tid;
		FloatingTypeGPU* hessian = forceAndHessian + 3;

		const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
		int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[vertexId * 2];
		int32_t numNeiTets = pTopology->vertexNeighborTets_infos[vertexId * 2 + 1];
		// i-th thread is in charge of i-th nei tet
		if (tid < numNeiTets)
		{
			int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + tid];
			int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
			CFloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * neiTetId;
			VEC12_SET(forceAndHessian, tetFAndHStorage);

		}
		else
		{
			VEC12_SET_CONSTANT(forceAndHessian, 0.f);
		}

		if (numNeiTets > NUM_THREADS_VERTEX_SWEEP)
			// if size larger than NUM_THREADS_VERTEX_SWEEP, we sequentially accumulate the extra elements to the elements of shared memory
		{
			const int extraSize = numNeiTets - NUM_THREADS_VERTEX_SWEEP;
			for (int i = 0; i < extraSize; i += NUM_THREADS_VERTEX_SWEEP)
			{
				int index = tid + NUM_THREADS_VERTEX_SWEEP + i;
				if (index < numNeiTets)
				{
					int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + index];
					int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
					CFloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * neiTetId;

					VEC12_ADD(forceAndHessian, tetFAndHStorage, forceAndHessian);
				}
			}
		}
		__syncthreads();

		for (unsigned int j = NUM_THREADS_VERTEX_SWEEP / 2; j > 0; j >>= 1)
		{
			if (tid < j) {
				float* forceAndHessianOther = forceAndHessianAll + (tid + j) * 12;
				VEC12_ADD(forceAndHessian, forceAndHessianOther, forceAndHessian);
			}

			__syncthreads();
		}

		// Unwarpped version

		//unwarppedReduction(forceAndHessianAll, tid);

		if (tid == 0)
		{
			accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysicsData->dtSqrReciprocal, forceAndHessian, hessian);
			accumulateCollisionForceAndHessian(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);
			accumulateBoundaryForceAndHessianGPU(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);

			if (CuMatrix::vec3NormSquare(hessian) > CMP_EPSILON2)
			{
				FloatingTypeGPU descentDirection[3];
				bool solverSuccess = CuMatrix::solve3x3_psd_stable(hessian, forceAndHessian, descentDirection);

				FloatingTypeGPU stepSize = solverSuccess ? pPhysicsData->stepSize : pPhysicsData->stepSizeGD;

				CuMatrix::vec3MulAddTo(descentDirection, stepSize, v);
			}

			if (pPhysicsData->useAccelerator
				&& acceleratorOmega != 1.0f)
			{
				// compute y_k+1
				FloatingTypeGPU* y_new = pTetMeshGPU->positionsNew + VERTEX_BUFFER_STRIDE * vertexId;
				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;
				if (pTetMeshGPU->activeCollisionMask[vertexId])
				{
					CuMatrix::vec3Set(y_new, v);
				}
				else
				{
					// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

					FloatingTypeGPU yDiff[3];
					CuMatrix::vec3Minus(v, y_prev, yDiff);
					CuMatrix::vec3Set(y_new, y_prev);
					CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, y_new);

				}
				// write y_k to y_k-1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}
	}
}



__global__ void GAIA::VBDSolveParallelGroup_allInOne_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega)
{
	__shared__ __builtin_align__(16) FloatingTypeGPU forceAndHessianAll[NUM_THREADS_VERTEX_SWEEP * 12];

	int32_t iV = blockIdx.x;
	int32_t meshId = parallelGroupsHead[iV * 2];
	int32_t vertexId = parallelGroupsHead[iV * 2 + 1];
	VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
	FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

	if (!pTetMeshGPU->vertexFixedMask[vertexId]
		&& pTetMeshGPU->activeForMaterialSolve
		//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
		)
	{
		int32_t tid = threadIdx.x;

		// assign material force & hessian to shared memory

		FloatingTypeGPU* forceAndHessian = forceAndHessianAll + 12 * tid;
		FloatingTypeGPU* hessian = forceAndHessian + 3;

		VEC12_SET_CONSTANT(forceAndHessian, 0.f);

		const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
		int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[vertexId * 2];
		int32_t numNeiTets = pTopology->vertexNeighborTets_infos[vertexId * 2 + 1];

		// i-th thread is in charge of i-th nei tet
		// evaluate the Hessian of the first 
		if (tid < numNeiTets)
		{
			int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + tid];
			int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + tid];
			int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;

			//CFloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * neiTetId;
			//VEC12_SET(forceAndHessian, tetFAndHStorage);
			// instead of using the precomputed hessians we compute them on the fly
			evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
				tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, forceAndHessian, hessian, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysicsData->dt);
		}

		if (numNeiTets > NUM_THREADS_VERTEX_SWEEP)
			// if size larger than NUM_THREADS_VERTEX_SWEEP, we sequentially accumulate the extra elements to the elements of shared memory
		{
			const int extraSize = numNeiTets - NUM_THREADS_VERTEX_SWEEP;
			for (int i = 0; i < extraSize; i += NUM_THREADS_VERTEX_SWEEP)
			{
				int index = tid + NUM_THREADS_VERTEX_SWEEP + i;
				if (index < numNeiTets)
				{
					int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + index];
					int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + index];
					int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;

					CFloatingTypeGPU tetRestVolume = pTetMeshGPU->tetRestVolume[neiTetId];
					// instead of using the precomputed hessians we compute them on the fly
					evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
						tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, forceAndHessian, hessian, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysicsData->dt);
				}
			}
		}
		__syncthreads();

		for (unsigned int j = NUM_THREADS_VERTEX_SWEEP / 2; j > 0; j >>= 1)
		{
			if (tid < j) {
				float* forceAndHessianOther = forceAndHessianAll + (tid + j) * 12;
				VEC12_ADD(forceAndHessian, forceAndHessianOther, forceAndHessian);
			}

			__syncthreads();
		}

		// Unwarpped version
		//unwarppedReduction(forceAndHessianAll, tid);

		if (tid == 0)
		{
			accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysicsData->dtSqrReciprocal, forceAndHessian, hessian);
			accumulateCollisionForceAndHessian(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);
			accumulateBoundaryForceAndHessianGPU(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);

			if (CuMatrix::vec3NormSquare(hessian) > CMP_EPSILON2)
			{
				FloatingTypeGPU descentDirection[3];
				bool solverSuccess = CuMatrix::solve3x3_psd_stable(hessian, forceAndHessian, descentDirection);

				FloatingTypeGPU stepSize = solverSuccess ? pPhysicsData->stepSize : pPhysicsData->stepSizeGD;

				CuMatrix::vec3MulAddTo(descentDirection, stepSize, v);
			}

			if (pPhysicsData->useAccelerator
				&& acceleratorOmega != 1.0f)
			{
				// compute y_k+1
				FloatingTypeGPU* y_new = pTetMeshGPU->positionsNew + VERTEX_BUFFER_STRIDE * vertexId;
				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;
				if (pTetMeshGPU->activeCollisionMask[vertexId])
				{
					CuMatrix::vec3Set(y_new, v);
				}
				else
				{
					// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

					FloatingTypeGPU yDiff[3];
					CuMatrix::vec3Minus(v, y_prev, yDiff);
					CuMatrix::vec3Set(y_new, y_prev);
					CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, y_new);

				}
				// write y_k to y_k-1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}
	}

	return ;
}

__global__ void GAIA::VBDSolveParallelGroup_vertexSweepAcceleratedGS_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize, FloatingTypeGPU acceleratorOmega)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < parallelGroupSize)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

		if (!pTetMeshGPU->vertexFixedMask[vertexId]
			&& pTetMeshGPU->activeForMaterialSolve
			//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
			)
		{
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);
			// compute y_k+1
			VBDStep_vertexSweep(pPhysicsData, pTetMeshGPU, vertexId);

			if (acceleratorOmega != 1.0f)
			{
				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;

				// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

				FloatingTypeGPU yDiff[3];
				CuMatrix::vec3Minus(v, y_prev, yDiff);
				CuMatrix::vec3Set(v, y_prev);
				CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, v);

				// write y_k to y_k+1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}

	}
	return;
}

__global__ void GAIA::GDSolveParallelGroup_vertexSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < parallelGroupSize)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		CFloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);

		if (!pTetMeshGPU->vertexFixedMask[vertexId]
			&& pTetMeshGPU->activeForMaterialSolve
			//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
			)
		{
			//if (pPhysicsData->useBlockJacobi)
			//{
			//	GDStep_vertexSweep_allInOne_blockJacobi(pPhysicsData, pTetMeshGPU, vertexId);
			//	//printf("BlockJacobi\n");
			//}
			//else {
			//	//printf("DiagonalJacobi\n");
			//}

			GDStep_vertexSweep_allInOne(pPhysicsData, pTetMeshGPU, vertexId);

		}

	}
	return;
}

__global__ void GAIA::GDSolveParallelGroup_BlockJacobi_allInOneSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, int32_t parallelGroupSize)
{
	__shared__ __builtin_align__(16) FloatingTypeGPU forceAndHessianAll[NUM_THREADS_VERTEX_SWEEP * 12];

	int32_t iV = blockIdx.x;
	int32_t meshId = parallelGroupsHead[iV * 2];
	int32_t vertexId = parallelGroupsHead[iV * 2 + 1];
	VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];

	if (!pTetMeshGPU->vertexFixedMask[vertexId]
		&& pTetMeshGPU->activeForMaterialSolve
		//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
		)
	{
		FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);
		FloatingTypeGPU* dx = pTetMeshGPU->positionsNew + 3 * vertexId;
		int32_t tid = threadIdx.x;

		// assign material force & hessian to shared memory

		FloatingTypeGPU* forceAndHessian = forceAndHessianAll + 12 * tid;
		FloatingTypeGPU* hessian = forceAndHessian + 3;

		VEC12_SET_CONSTANT(forceAndHessian, 0.f);

		const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
		int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[vertexId * 2];
		int32_t numNeiTets = pTopology->vertexNeighborTets_infos[vertexId * 2 + 1];

		// i-th thread is in charge of i-th nei tet
		// evaluate the Hessian of the first 
		if (tid < numNeiTets)
		{
			int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + tid];
			int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + tid];
			int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;

			//CFloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * neiTetId;
			//VEC12_SET(forceAndHessian, tetFAndHStorage);
			// instead of using the precomputed hessians we compute them on the fly
			evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
				tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, forceAndHessian, hessian, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysicsData->dt);
		}

		if (numNeiTets > NUM_THREADS_VERTEX_SWEEP)
			// if size larger than NUM_THREADS_VERTEX_SWEEP, we sequentially accumulate the extra elements to the elements of shared memory
		{
			const int extraSize = numNeiTets - NUM_THREADS_VERTEX_SWEEP;
			for (int i = 0; i < extraSize; i += NUM_THREADS_VERTEX_SWEEP)
			{
				int index = tid + NUM_THREADS_VERTEX_SWEEP + i;
				if (index < numNeiTets)
				{
					int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + index];
					int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + index];
					int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;

					CFloatingTypeGPU tetRestVolume = pTetMeshGPU->tetRestVolume[neiTetId];
					// instead of using the precomputed hessians we compute them on the fly
					evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
						tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, forceAndHessian, hessian, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysicsData->dt);
				}
			}
		}
		__syncthreads();

		for (unsigned int j = NUM_THREADS_VERTEX_SWEEP / 2; j > 0; j >>= 1)
		{
			if (tid < j) {
				float* forceAndHessianOther = forceAndHessianAll + (tid + j) * 12;
				VEC12_ADD(forceAndHessian, forceAndHessianOther, forceAndHessian);
			}

			__syncthreads();
		}

		// Unwarpped version
		//unwarppedReduction(forceAndHessianAll, tid);

		if (tid == 0)
		{
			accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysicsData->dtSqrReciprocal, forceAndHessian, hessian);
			accumulateCollisionForceAndHessian(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);
			accumulateBoundaryForceAndHessianGPU(pPhysicsData, pTetMeshGPU, vertexId, forceAndHessian, hessian);

			if (CuMatrix::vec3NormSquare(hessian) > CMP_EPSILON2)
			{
				FloatingTypeGPU descentDirection[3];
				bool solverSuccess = CuMatrix::solve3x3_psd_stable(hessian, forceAndHessian, descentDirection);

#ifdef PRINT_DBG_INFO
				printf("descentDirection: ");
				CuMatrix::printFloatVec(descentDirection, 3);
#endif
#ifdef GPU_JACOBI_DX
				if (solverSuccess)
				{
					CuMatrix::vec3Set(dx, descentDirection);
				}
				else
				{
					CuMatrix::vec3Set(dx, 0.f);
				}
#endif // GPU_JACOBI_DX

			}
			else
			{
				CuMatrix::vec3Set(dx, 0.f);

			}
		}
		else
		{
			CuMatrix::vec3Set(dx, 0.f);

		}
	}

	return ;
}

__global__ void GAIA::GDSolveParallelGroup_updatePositionSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, CFloatingTypeGPU stepSize, CFloatingTypeGPU acceleratorOmega, int32_t* parallelGroupsHead, int32_t parallelGroupSize)
{
	int32_t iVert = blockIdx.x * blockDim.x + threadIdx.x;
	if (iVert < parallelGroupSize)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		FloatingTypeGPU* v = pTetMeshGPU->getVert(vertexId);
		CFloatingTypeGPU* dx = pTetMeshGPU->positionsNew + VERTEX_BUFFER_STRIDE * vertexId;

		if (!pTetMeshGPU->vertexFixedMask[vertexId]
			&& pTetMeshGPU->activeForMaterialSolve
			//&& v[GRAVITY_DIM] < pPhysicsData->collisionOffHeight
			)
		{
			FloatingTypeGPU y_k[3];
			// record y_k
			CuMatrix::vec3Set(y_k, v);
			CuMatrix::vec3MulAddTo(dx, stepSize, v);

			if (acceleratorOmega != 1.0f)
			{
				FloatingTypeGPU* y_prev = pTetMeshGPU->positionsPrevIter + VERTEX_BUFFER_STRIDE * vertexId;

				// y_k+1_accelerated = (y_k+1 - y_k-1) * omega +  y_k-1

				FloatingTypeGPU yDiff[3];
				CuMatrix::vec3Minus(v, y_prev, yDiff);
				CuMatrix::vec3Set(v, y_prev);
				CuMatrix::vec3MulAddTo(yDiff, acceleratorOmega, v);

				// write y_k to y_k-1
				CuMatrix::vec3Set(y_prev, y_k);
			}
		}

	}
}

__global__ void GAIA::VBDSolveParallelGroup_kernel_serial(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	int32_t parallelGroupSize)
{
	for (int32_t iVert = 0; iVert < parallelGroupSize; iVert++)
	{
		int32_t meshId = parallelGroupsHead[iVert * 2];
		int32_t vertexId = parallelGroupsHead[iVert * 2 + 1];
		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
#ifdef PRINT_DBG_INFO
		printf("-------------------------------------\nProcessing vertex %d of mesh %d\n", vertexId, meshId);
#endif
		if (!pTetMeshGPU->vertexFixedMask[vertexId])
		{
			VBDStep_allInOne(pPhysicsData, pTetMeshGPU, vertexId);
		}
	}
}

__global__ void GAIA::evaluateElasticEnergy_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize, FloatingTypeGPU* tetElasticEnergy)
{
	int32_t iTet = blockIdx.x * blockDim.x + threadIdx.x;
	if (iTet < tetParallelGroupSize)
	{
		int32_t meshId = tetParallelGroupsHead[iTet * 2];
		int32_t tetId = tetParallelGroupsHead[iTet * 2 + 1];

		VBDTetMeshNeoHookeanGPU* pTetMeshGPU = (VBDTetMeshNeoHookeanGPU*)pPhysicsData->tetMeshes[meshId];
		int32_t* tetVIds = pTetMeshGPU->pTopology->tetVIds + 4 * tetId;

		tetElasticEnergy[iTet] = evaluateNeoHookeanEnergy(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[tetId], tetId, 
			tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->DmInvs);
	}
}



GPU_CPU_INLINE_FUNC void GAIA::caculateMaterialForceAndHessian_NeoHookean_oneTet(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t tetId, int32_t vertexOrderInTet,
	FloatingTypeGPU* force, FloatingTypeGPU* h, CFloatingTypeGPU dt)
{
	const int32_t* tetVIds = pTetMeshGPU->pTopology->tetVIds + 4 * tetId;
	evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[tetId], tetId, vertexOrderInTet,
		tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, force, h + 3, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, dt);
	return;
}

GPU_CPU_INLINE_FUNC void GAIA::updateVertexVelocity(VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU dt)
{
	CFloatingTypeGPU* pos = pTetMeshGPU->getVert(vertexId);
	CFloatingTypeGPU* posPrev = pTetMeshGPU->vertPrevPos + VERTEX_BUFFER_STRIDE * vertexId;
	FloatingTypeGPU* vel = pTetMeshGPU->velocity + VERTEX_BUFFER_STRIDE * vertexId;

	CuMatrix::vec3Minus(pos, posPrev, vel);
	CuMatrix::vec3Mul(vel, 1.f / dt, vel);

	return;
}

GPU_CPU_INLINE_FUNC void GAIA::applyBoundaryFriction(VBDPhysicsDataGPU* pPhysics, VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId,
	CFloatingTypeGPU dt, CFloatingTypeGPU vertInvMass)
{
	CFloatingTypeGPU boundaryCollisionStiffness = pPhysics->boundaryCollisionStiffness;
	CFloatingTypeGPU frictionRatio = pPhysics->boundaryFrictionDynamic;
	CFloatingTypeGPU* pos = pTetMeshGPU->getVert(vertexId);
	FloatingTypeGPU* vel = pTetMeshGPU->velocity + VERTEX_BUFFER_STRIDE * vertexId;

	for (int32_t iDim = 0; iDim < 3; iDim++)
	{
		FloatingTypeGPU contactNormal[3] = { 0.f, 0.f, 0.f };

		CFloatingTypeGPU lowerBound = pPhysics->worldBounds[iDim];
		CFloatingTypeGPU upperBound = pPhysics->worldBounds[iDim + 3];

		if (pos[iDim] < lowerBound)
		{
			CFloatingTypeGPU penetrationDepth = lowerBound - pos[iDim];

			contactNormal[iDim] = 1;
			applyFrictionalVelocityDampingBoundary(vel, contactNormal, vertInvMass, frictionRatio, penetrationDepth * boundaryCollisionStiffness, dt);
		}
		else if (pos[iDim] > upperBound)
		{
			CFloatingTypeGPU penetrationDepth = pos[iDim] - upperBound;

			contactNormal[iDim] = -1;
			applyFrictionalVelocityDampingBoundary(vel, contactNormal, vertInvMass, frictionRatio, penetrationDepth * boundaryCollisionStiffness, dt);
		}
	}
	return;
}

GPU_CPU_INLINE_FUNC void GAIA::applyCollisionFriction(VBDPhysicsDataGPU* pPhysics, VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU dt, CFloatingTypeGPU vertInvMass)
{
	FloatingTypeGPU* vel = pTetMeshGPU->velocity + VERTEX_BUFFER_STRIDE * vertexId;
	//if (pTetMeshGPU->activeCollisionMask[vertexId])
	//{
	CFloatingTypeGPU frictionRatio = pTetMeshGPU->frictionDynamic;
	FloatingTypeGPU contactNormal[3];

	CollisionDataGPU& collisionData = pTetMeshGPU->getCollisionData(vertexId);

	if (collisionData.activeColliding)
	{
		// v side of the v-f collisions
		CFloatingTypeGPU* tangentialRelativeVelDirection = collisionData.collisionForceAndHessian;
		CFloatingTypeGPU* tangentialRelativeVelNorm = collisionData.collisionForceAndHessian + 3;
		CFloatingTypeGPU contactForceNorm = *(collisionData.collisionForceAndHessian + 4);
		if (contactForceNorm)
		{
			FloatingTypeGPU frictionForceDirection[3];
			CuMatrix::vec3Mul(tangentialRelativeVelDirection, -1.f, frictionForceDirection);
			applyFrictionalVelocityDamping(vel, frictionForceDirection, *tangentialRelativeVelNorm,
				contactForceNorm, vertInvMass, frictionRatio, dt);
		}
	}

	// f side of the v-f collision
	for (size_t iCollision = 0; iCollision < collisionData.numCollisionRelations; iCollision++)
	{
		CollisionRelationGPU& colRelation = collisionData.collisionRelations[iCollision];
		VBDBaseTetMeshGPU* pTetMeshGPU_other = pPhysics->tetMeshes[colRelation.meshId];

		CollisionDataGPU& colResultOther = pTetMeshGPU_other->getCollisionData(colRelation.collisionPrimitiveId);

		/*	printf("collision handling for vertex: %d, collisionPrimitiveId: d, collisionOrder: %d\n",
				iV, collisionRelation.collisionPrimitiveId, collisionRelation.collisionVertexOrder);
			printf("    collision info:\n");
			printf("    collide with (%d, %d, %d) from mesh %d:\n",
				collisionDataOther.closestSurfaceFaceVIds[0], collisionDataOther.closestSurfaceFaceVIds[1], collisionDataOther.closestSurfaceFaceVIds[2],
				collisionDataOther.intersectedTMeshId);*/
		CFloatingTypeGPU* tangentialRelativeVelDirection = colResultOther.collisionForceAndHessian;
		CFloatingTypeGPU* tangentialRelativeVelNorm = colResultOther.collisionForceAndHessian + 3;
		CFloatingTypeGPU b = colResultOther.closestSurfacePtBarycentrics[colRelation.collisionVertexOrder];
		CFloatingTypeGPU contactForceNorm = *(colResultOther.collisionForceAndHessian + 4) * b;
		if (contactForceNorm)
		{
			applyFrictionalVelocityDamping(vel, tangentialRelativeVelDirection, *tangentialRelativeVelNorm,
				contactForceNorm, vertInvMass, frictionRatio, dt);
		}

	}
	//}
	return;
}

GPU_CPU_INLINE_FUNC void GAIA::applyFrictionalVelocityDamping(FloatingTypeGPU* velocity, CFloatingTypeGPU* frictionForceDirection,
	CFloatingTypeGPU tangentialVelNorm, CFloatingTypeGPU contactForceNorm, CFloatingTypeGPU vertInvMass, CFloatingTypeGPU frictionRatio, CFloatingTypeGPU dt)
{
	//FloatingTypeGPU tangentialVel[3];
	//FloatingTypeGPU tangentialVelMag = CuMatrix::vec3DotProduct(velocity, tangentialForceDirection);
	//CuMatrix::vec3Mul(tangentialForceDirection, tangentialVelMag, tangentialVel);

	//FloatingTypeGPU contactVel[3];
	//CuMatrix::vec3Minus(velocity, tangentialVel, contactVel);
	if (abs(tangentialVelNorm) > 1e-6f)
	{
		FloatingTypeGPU velDamping = contactForceNorm * vertInvMass * dt * frictionRatio;
		velDamping = velDamping <= 0.5f * tangentialVelNorm ? velDamping : 0.5f * tangentialVelNorm;

		CuMatrix::vec3MulAddTo(frictionForceDirection, velDamping, velocity);
	}

}

GPU_CPU_INLINE_FUNC void GAIA::applyFrictionalVelocityDampingBoundary(FloatingTypeGPU* velocity, CFloatingTypeGPU* contactNormal, CFloatingTypeGPU vertInvMass,
	CFloatingTypeGPU frictionRatio, CFloatingTypeGPU contactForce, CFloatingTypeGPU dt)
{
	FloatingTypeGPU contactVel[3];
	FloatingTypeGPU contactVelMag = CuMatrix::vec3DotProduct(velocity, contactNormal);
	CuMatrix::vec3Mul(contactNormal, contactVelMag, contactVel);

	FloatingTypeGPU orthogonalVel[3];
	CuMatrix::vec3Minus(velocity, contactVel, orthogonalVel);
	CFloatingTypeGPU orthogonalVelNorm = CuMatrix::vec3Norm(orthogonalVel);
	if (orthogonalVelNorm > 1e-6f)
	{
		CFloatingTypeGPU velDamping = contactForce * vertInvMass * dt * frictionRatio;

		if (velDamping >= orthogonalVelNorm)
		{
			CuMatrix::vec3Set(orthogonalVel, 0.f);
		}
		else {
			CFloatingTypeGPU orthogonalVelNormDamped = orthogonalVelNorm - velDamping;
			// orthogonalVel = orthogonalVel * (orthogonalVelNormDamped / orthogonalVelNorm);
			CuMatrix::vec3Mul(orthogonalVel, orthogonalVelNormDamped / orthogonalVelNorm, orthogonalVel);
		}
	}

	CuMatrix::vec3Add(orthogonalVel, contactVel, velocity);

	return;
}

GPU_CPU_INLINE_FUNC void GAIA::dampVelocity(FloatingTypeGPU* velocity, CFloatingTypeGPU maxVelocityMagnitude, CFloatingTypeGPU exponentialVelDamping, CFloatingTypeGPU constantVelDamping)
{
	CFloatingTypeGPU vMag = CuMatrix::vec3Norm(velocity);

	if (maxVelocityMagnitude > 0 && vMag > maxVelocityMagnitude)
	{
		CuMatrix::vec3Mul(velocity, maxVelocityMagnitude / vMag, velocity);
	}
	else if (vMag > 1e-6f)
	{
		FloatingTypeGPU vMagNew = vMag * exponentialVelDamping - constantVelDamping;
		vMagNew = vMagNew > 1e-6f ? vMagNew : 0.f;
		CuMatrix::vec3Mul(velocity, vMagNew / vMag, velocity);
	}

	return;
}

__device__ __host__ void GAIA::updateCollisionInfoGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t meshId, int32_t vertexId)
{
	CollisionDataGPU& collisionData = pPhysicsData->tetMeshes[meshId]->getCollisionData(vertexId);

	/*printf("Sweeping collision %d from mesh %d\n", vertexId, meshId);
	printVBDPhysicsDataGPU(&collisionData);*/

	const VBDBaseTetMeshGPU* pTetMeshGPU = pPhysicsData->tetMeshes[meshId];
	CFloatingTypeGPU* p = pTetMeshGPU->getVert(vertexId);

	const VBDBaseTetMeshGPU* pTetMeshGPU_colliding = pPhysicsData->tetMeshes[collisionData.intersectedTMeshId];

	// update the closest point
	CFloatingTypeGPU* a = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[0]);
	CFloatingTypeGPU* b = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[1]);
	CFloatingTypeGPU* c = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[2]);

	FloatingTypeGPU* closestSurfacePtBarycentrics = collisionData.closestSurfacePtBarycentrics;
	FloatingTypeGPU closestSurfacePt[3];

	//ClosestPointOnTriangleTypeGPU closestPtType =
	//	closestPointTriangle(p, a, b, c, collisionData.closestSurfacePtBarycentrics, collisionData.closestSurfacePt);

	ClosestPointOnTriangleTypeGPU closestPtType =
		closestPointTriangle(p, a, b, c, closestSurfacePtBarycentrics, closestSurfacePt);

	// compute contact normal
	// face interior
	//FloatingTypeGPU* normal = collisionData.closestPointNormal;
	FloatingTypeGPU normal[3];
	FloatingTypeGPU faceNormal[3];
	CuMatrix::faceNormal(a, b, c, faceNormal);
	FloatingTypeGPU diff[3];
	CuMatrix::vec3Minus(closestSurfacePt, p, diff);

	if (closestPtType == GAIA::ClosestPointOnTriangleTypeGPU::AtInterior)
	{
		CuMatrix::vec3Set(normal, faceNormal);
	}
	// all other cases
	else if (closestPtType != GAIA::ClosestPointOnTriangleTypeGPU::NotFound)
	{
		CuMatrix::vec3Mul(diff, 1.f / CuMatrix::vec3Norm(diff), normal);

		if (CuMatrix::vec3DotProduct(normal, faceNormal) < 0)
		{
			CuMatrix::vec3Mul(normal, -1.f, normal);
		}
	}
	CFloatingTypeGPU penetrationDepth = CuMatrix::vec3DotProduct(diff, normal);

	// compute the force and hessian
	for (size_t i = 0; i < 12; i++)
	{
		collisionData.collisionForceAndHessian[i] = 0.f;
	}
	if (penetrationDepth > 0.f)
	{
		// collision force
		CFloatingTypeGPU k = pPhysicsData->collisionStiffness;
		CFloatingTypeGPU lambda = penetrationDepth * k;
		FloatingTypeGPU* collisionForce = collisionData.collisionForceAndHessian;
		FloatingTypeGPU* collisionHessian = collisionData.collisionForceAndHessian + 3;

		CuMatrix::vec3MulAddTo(normal, lambda, collisionForce);
		CuMatrix::vec3OuterProduct(normal, normal, collisionHessian);
		CuMatrix::vec9MulAddTo(collisionHessian, k, collisionHessian);

#ifdef USE_IPC_FRICTION
		// friction
		FloatingTypeGPU dx3[3];
		CuMatrix::vec3Minus(p, pTetMeshGPU->getVertPrevPos(vertexId), dx3);
		FloatingTypeGPU dx0[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[0]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[0]), dx0);

		//pIntersectedTM->vertex(t0) - pIntersectedTM->vertexPrevPos(t0);
		FloatingTypeGPU dx1[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[1]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[1]), dx1);
		FloatingTypeGPU dx2[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[2]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[2]), dx2);

		CFloatingTypeGPU bary0 = closestSurfacePtBarycentrics[0];
		CFloatingTypeGPU bary1 = closestSurfacePtBarycentrics[1];
		CFloatingTypeGPU bary2 = closestSurfacePtBarycentrics[2];
		CFloatingTypeGPU dx[3] = {
			dx3[0] - (bary0 * dx0[0] + bary1 * dx1[0] + bary2 * dx2[0]),
			dx3[1] - (bary0 * dx0[1] + bary1 * dx1[1] + bary2 * dx2[1]),
			dx3[2] - (bary0 * dx0[2] + bary1 * dx1[2] + bary2 * dx2[2])
		};
		FloatingTypeGPU e0[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[1]),
			pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[0]),
			e0);
		CuMatrix::vec3Normalize(e0);
		//e0 = (pIntersectedTM->vertex(t1) - pIntersectedTM->vertex(t0)).normalized();

		FloatingTypeGPU e1[3];
		CuMatrix::vec3CrossProduct(e0, normal, e1);
		CuMatrix::vec3Normalize(e1);

		FloatingTypeGPU T[6] = {
			e0[0],
			e0[1],
			e0[2],
			e1[0],
			e1[1],
			e1[2]
		};
		//T.col(0) = e0;
		//T.col(1) = (e0.cross(n)).normalized();

		CFloatingTypeGPU u[2] = {
		T[0] * dx[0] + T[1] * dx[1] + T[2] * dx[2],
		T[3] * dx[0] + T[4] * dx[1] + T[5] * dx[2],
		};
		// Vec2 u = T.transpose() * positionsNew;

		// average of the two friction coefficients
		CFloatingTypeGPU mu = (pTetMeshGPU->frictionDynamic + pTetMeshGPU_colliding->frictionDynamic) * 0.5f;
		CFloatingTypeGPU epsV = (pTetMeshGPU->frictionEpsV + pTetMeshGPU_colliding->frictionEpsV) * 0.5f;
		CFloatingTypeGPU dt = pPhysicsData->dt;
		CFloatingTypeGPU epsU = epsV * dt;

		accumulateVertexFrictionGPU(mu, lambda, T, u, epsU, collisionForce, collisionHessian);
#endif // USE_IPC_FRICTION
	}

	/*printf("normal %f %f %f\n", normal[0], normal[1], normal[2]);
	printf("closestSurfacePt %f %f %f\n",
		collisionData.closestSurfacePt[0], collisionData.closestSurfacePt[1], collisionData.closestSurfacePt[2]);*/


		//switch (closestPtType)
		//{
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtA:
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtB:
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtC:
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtAB:
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtBC:
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtAC:
		//	FloatingTypeGPU* normal = collisionData.closestPointNormal;
		//	CuMatrix::vec3Minus(collisionData.closestSurfacePt, p, normal);
		//	CuMatrix::vec3Normalize(normal);
		//	break;
		//case GAIA::ClosestPointOnTriangleTypeGPU::AtInterior:
		//	pIntersectedTM->computeFaceNormal(closestFaceId, normal);
		//	break;
		//case GAIA::ClosestPointOnTriangleTypeGPU::NotFound:
		//	return;
		//	break;
		//default:
		//	return;
		//	break;
		//}
	return;
}

__device__ __host__ void GAIA::computeRelativeVelocityGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t meshId, int32_t vertexId)
{
	CollisionDataGPU& collisionData = pPhysicsData->tetMeshes[meshId]->getCollisionData(vertexId);
	if (!collisionData.activeColliding)
	{
		return;
	}
	/*printf("Sweeping collision %d from mesh %d\n", vertexId, meshId);
	printVBDPhysicsDataGPU(&collisionData);*/

	const VBDBaseTetMeshGPU* pTetMeshGPU = pPhysicsData->tetMeshes[meshId];
	CFloatingTypeGPU* p = pTetMeshGPU->getVert(vertexId);

	const VBDBaseTetMeshGPU* pTetMeshGPU_colliding = pPhysicsData->tetMeshes[collisionData.intersectedTMeshId];

	// update the closest point
	CFloatingTypeGPU* a = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[0]);
	CFloatingTypeGPU* b = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[1]);
	CFloatingTypeGPU* c = pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[2]);

	FloatingTypeGPU* closestSurfacePtBarycentrics = collisionData.closestSurfacePtBarycentrics;
	FloatingTypeGPU closestSurfacePt[3];

	//ClosestPointOnTriangleTypeGPU closestPtType =
	//	closestPointTriangle(p, a, b, c, collisionData.closestSurfacePtBarycentrics, collisionData.closestSurfacePt);

	ClosestPointOnTriangleTypeGPU closestPtType =
		closestPointTriangle(p, a, b, c, closestSurfacePtBarycentrics, closestSurfacePt);

	// compute contact normal
	// face interior
	//FloatingTypeGPU* normal = collisionData.closestPointNormal;
	FloatingTypeGPU normal[3];
	FloatingTypeGPU faceNormal[3];
	CuMatrix::faceNormal(a, b, c, faceNormal);
	FloatingTypeGPU diff[3];
	CuMatrix::vec3Minus(closestSurfacePt, p, diff);
	if (closestPtType == GAIA::ClosestPointOnTriangleTypeGPU::AtInterior)
	{
		CuMatrix::vec3Set(normal, faceNormal);
	}
	// all other cases
	else if (closestPtType != GAIA::ClosestPointOnTriangleTypeGPU::NotFound)
	{
		CuMatrix::vec3Mul(diff, 1.f / CuMatrix::vec3Norm(diff), normal);

		if (CuMatrix::vec3DotProduct(normal, faceNormal) < 0)
		{
			CuMatrix::vec3Mul(normal, -1.f, normal);
		}
	}
	CFloatingTypeGPU penetrationDepth = CuMatrix::vec3DotProduct(diff, normal);

	// compute the force and hessian
	for (size_t i = 0; i < 12; i++)
	{
		collisionData.collisionForceAndHessian[i] = 0.f;
	}

	FloatingTypeGPU* tangentialRelativeVelDirection = collisionData.collisionForceAndHessian;
	FloatingTypeGPU* tangentialRelativeVelNorm = collisionData.collisionForceAndHessian + 3;
	FloatingTypeGPU* contactForceNorm = collisionData.collisionForceAndHessian + 4;

	if (penetrationDepth > 0.f)
	{
		// collision force
		CFloatingTypeGPU k = pPhysicsData->collisionStiffness;
		CFloatingTypeGPU lambda = penetrationDepth * k;

		// friction
		FloatingTypeGPU dx3[3];
		CuMatrix::vec3Minus(p, pTetMeshGPU->getVertPrevPos(vertexId), dx3);
		FloatingTypeGPU dx0[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[0]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[0]), dx0);
		FloatingTypeGPU dx1[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[1]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[1]), dx1);
		FloatingTypeGPU dx2[3];
		CuMatrix::vec3Minus(pTetMeshGPU_colliding->getVert(collisionData.closestSurfaceFaceVIds[2]),
			pTetMeshGPU_colliding->getVertPrevPos(collisionData.closestSurfaceFaceVIds[2]), dx2);
		CFloatingTypeGPU bary0 = closestSurfacePtBarycentrics[0];
		CFloatingTypeGPU bary1 = closestSurfacePtBarycentrics[1];
		CFloatingTypeGPU bary2 = closestSurfacePtBarycentrics[2];
		CFloatingTypeGPU dx[3] = {
			dx3[0] - (bary0 * dx0[0] + bary1 * dx1[0] + bary2 * dx2[0]),
			dx3[1] - (bary0 * dx0[1] + bary1 * dx1[1] + bary2 * dx2[1]),
			dx3[2] - (bary0 * dx0[2] + bary1 * dx1[2] + bary2 * dx2[2])
		};
		CFloatingTypeGPU dt = pPhysicsData->dt;

		CFloatingTypeGPU relativeVelocity[3] = {
			dx[0] / dt,
			dx[1] / dt,
			dx[2] / dt,
		};

		FloatingTypeGPU orthogonalVel[3];
		FloatingTypeGPU contactVel[3];

		FloatingTypeGPU contactVelMag = CuMatrix::vec3DotProduct(relativeVelocity, normal);
		CuMatrix::vec3Mul(normal, contactVelMag, contactVel);
		CuMatrix::vec3Minus(relativeVelocity, contactVel, orthogonalVel);
		CFloatingTypeGPU orthogonalVelNorm = CuMatrix::vec3Norm(orthogonalVel);
		CuMatrix::vec3Mul(orthogonalVel, 1.f / orthogonalVelNorm, orthogonalVel);

		CuMatrix::vec3Set(tangentialRelativeVelDirection, orthogonalVel);
		*tangentialRelativeVelNorm = orthogonalVelNorm;
		*contactForceNorm = lambda;
	}
	else
	{
		*tangentialRelativeVelNorm = 0.f;
		*contactForceNorm = 0.f;
	}

	return;
}

__device__ __host__  void GAIA::accumulateMaterialForceAndHessianForVertex_NeoHookean(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t iV,
	FloatingTypeGPU* force, FloatingTypeGPU* h, CFloatingTypeGPU dt)
{
	const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
	int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[iV * 2];
	int32_t numNeiTets = pTopology->vertexNeighborTets_infos[iV * 2 + 1];

	for (int32_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++)
	{
		int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + iNeiTet];
		int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
		int32_t vOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + iNeiTet];
		evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vOrderInTet,
			tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, force, h, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, dt);
	}
}

GPU_CPU_INLINE_FUNC void GAIA::accumulateMaterialForceAndHessianForVertex_NeoHookean_preCompuated(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t iV, FloatingTypeGPU* force,
	FloatingTypeGPU* h)
{
	const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
	int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[iV * 2];
	int32_t numNeiTets = pTopology->vertexNeighborTets_infos[iV * 2 + 1];

	for (int32_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++)
	{
		int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + iNeiTet];
		int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
		CFloatingTypeGPU* tetFAndHStorage = pTetMeshGPU->tetForceAndHessians + 12 * neiTetId;

		CuMatrix::vec3Add(tetFAndHStorage, force, force);
		CuMatrix::vec9Add(tetFAndHStorage + 3, h, h);

	}
	return;
}

__device__ __host__ void GAIA::accumulateCollisionForceAndHessian(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
	int32_t iV, FloatingTypeGPU* force, FloatingTypeGPU* h)
{
	const CollisionDataGPU& collisionData = pTetMeshGPU->getCollisionData(iV);

	bool activeColliding = false;
	if (collisionData.activeColliding)
	{
		//printf("collision handling for vertex: %d, collisionOrder: %d\n", iV, 3);
		//printf("    collision info:\n");
		//printf("    collide with (%d, %d, %d) from mesh %d:\n",
		//	collisionData.closestSurfaceFaceVIds[0], collisionData.closestSurfaceFaceVIds[1], collisionData.closestSurfaceFaceVIds[2],
		//	collisionData.intersectedTMeshId);

		accumlateCollisionForceAndHessianPerCollision(pPhysicsData, pTetMeshGPU, collisionData, iV, 3, force, h);
		activeColliding = true;
	}

	if (collisionData.numCollisionRelations)
	{
		for (int32_t iCollision = 0; iCollision < collisionData.numCollisionRelations; iCollision++)
		{
			const CollisionRelationGPU& collisionRelation = collisionData.collisionRelations[iCollision];
			const VBDBaseTetMeshGPU* pTetMeshGPUOther = pPhysicsData->tetMeshes[collisionRelation.meshId];
			const CollisionDataGPU& collisionDataOther = pTetMeshGPUOther->getCollisionData(collisionRelation.collisionPrimitiveId);

			/*	printf("collision handling for vertex: %d, collisionPrimitiveId: d, collisionOrder: %d\n",
					iV, collisionRelation.collisionPrimitiveId, collisionRelation.collisionVertexOrder);
				printf("    collision info:\n");
				printf("    collide with (%d, %d, %d) from mesh %d:\n",
					collisionDataOther.closestSurfaceFaceVIds[0], collisionDataOther.closestSurfaceFaceVIds[1], collisionDataOther.closestSurfaceFaceVIds[2],
					collisionDataOther.intersectedTMeshId);*/

			accumlateCollisionForceAndHessianPerCollision(pPhysicsData, pTetMeshGPUOther, collisionDataOther,
				collisionRelation.collisionPrimitiveId, collisionRelation.collisionVertexOrder, force, h);
			activeColliding = true;
		}
	}

	pTetMeshGPU->activeCollisionMask[iV] = activeColliding;

	return;
}

GPU_CPU_INLINE_FUNC CFloatingTypeGPU GAIA::evaluateNeoHookeanEnergy(FloatingTypeGPU miu, FloatingTypeGPU lmbd, FloatingTypeGPU A, int32_t tetId,
	const int32_t* tetVIds, FloatingTypeGPU* vs, CFloatingTypeGPU* DmInvs)
{
	CFloatingTypeGPU* v0 = vs + 3 * tetVIds[0];
	CFloatingTypeGPU* v1 = vs + 3 * tetVIds[1];
	CFloatingTypeGPU* v2 = vs + 3 * tetVIds[2];
	CFloatingTypeGPU* v3 = vs + 3 * tetVIds[3];

	FloatingTypeGPU Ds[9];
	FloatingTypeGPU F[9];

	CuMatrix::vec3Minus(v1, v0, Ds);
	CuMatrix::vec3Minus(v2, v0, Ds + 3);
	CuMatrix::vec3Minus(v3, v0, Ds + 6);
	CFloatingTypeGPU* DmInv = DmInvs + tetId * 9;
	//printf("Ds:\n");
	//CuMatrix::printMat3(Ds);

	CuMatrix::mat3MatProduct(Ds, DmInv, F);
	//printf("F:\n");
	//CuMatrix::printMat3(F);

	CFloatingTypeGPU a = 1 + miu / lmbd;

	CFloatingTypeGPU Phi_D = 0.5f * (CuMatrix::mat3FNormSquare(F) - 3.f);
	CFloatingTypeGPU detF = CuMatrix::mat3Determinant(F);
	CFloatingTypeGPU Phi_H = 0.5f * SQR(detF - a);

	FloatingTypeGPU E = A * (miu * Phi_D + lmbd * Phi_H);

	CFloatingTypeGPU restEnergy = 0.5 * A * lmbd * SQR(1.f - a);

	E -= restEnergy;

	return E;
}

GPU_CPU_INLINE_FUNC  FloatingTypeGPU GAIA::computeCollisionForcePerCollision(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
	const CollisionDataGPU& collisionResult, int32_t collisionResultId, int32_t collisionVertexOrder, FloatingTypeGPU* contactNormal)
{
	//int collidedMeshID = collisionResult.intersectedTMeshId;
	//const VBDBaseTetMeshGPU* pIntersectedTM = pPhysicsData->tetMeshes[collidedMeshID];
	////CFloatingTypeGPU* p = pTetMeshGPU->getVert(collisionResultId);
	////CFloatingTypeGPU k = pPhysicsData->collisionStiffness;
	//CFloatingTypeGPU collisionRadius = pPhysicsData->collisionAirDistance;
	////CFloatingTypeGPU* n = collisionResult.closestPointNormal;
	////CFloatingTypeGPU* c = collisionResult.closestSurfacePt;
	////FloatingTypeGPU diff[3];
	////CuMatrix::vec3Minus(c, p, diff);

	//CFloatingTypeGPU penetrationDepth = CuMatrix::vec3DotProduct(diff, n);
	//FloatingTypeGPU collisionHessian[9];
	//if (penetrationDepth > 0)
	//{
	//	FloatingTypeGPU b;
	//	switch (collisionVertexOrder)
	//	{
	//	case 0:
	//		b = -collisionResult.closestSurfacePtBarycentrics[0];
	//		break;
	//	case 1:
	//		b = -collisionResult.closestSurfacePtBarycentrics[1];
	//		break;
	//	case 2:
	//		b = -collisionResult.closestSurfacePtBarycentrics[2];
	//		break;
	//	case 3:
	//		b = 1.0f;
	//		// CuMatrix::vec3Mul(diff, -1.f, diff);
	//		break;
	//	default:
	//		break;
	//	}

	//	if (b != 0)
	//	{
	//		// simplified point-plane
	//		CuMatrix::vec3Set(contactNormal, n);
	//		return abs(k * b * penetrationDepth);
	//	}
	//	else
	//	{
	//		return 0.f;
	//	}
	//}
	return 0.f;
}

GPU_CPU_INLINE_FUNC void GAIA::accumlateCollisionForceAndHessianPerCollision(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
	const CollisionDataGPU& collisionResult, int32_t collisionResultId,
	int32_t collisionVertexOrder, FloatingTypeGPU* force, FloatingTypeGPU* hessian)
{
	int collidedMeshID = collisionResult.intersectedTMeshId;

	FloatingTypeGPU b;
	switch (collisionVertexOrder)
	{
	case 0:
		b = -collisionResult.closestSurfacePtBarycentrics[0];
		break;
	case 1:
		b = -collisionResult.closestSurfacePtBarycentrics[1];
		break;
	case 2:
		b = -collisionResult.closestSurfacePtBarycentrics[2];
		break;
	case 3:
		b = 1.0f;
		// CuMatrix::vec3Mul(diff, -1.f, diff);
		break;
	default:
		break;
	}

	// collision force and hessian is precomputed in updateCollisionInfoGPU
	CuMatrix::vec3MulAddTo(collisionResult.collisionForceAndHessian, b, force);
	CuMatrix::vec9MulAddTo(collisionResult.collisionForceAndHessian + 3, b * b, hessian);

	// simplified point-plane
	//CuMatrix::vec3MulAddTo(n, k * b * penetrationDepth, force);
	//CuMatrix::vec3OuterProduct(n, n, collisionHessian);
	//CuMatrix::vec9MulAddTo(collisionHessian, k * b * b, hessian);

	//printf("after accumulatCollisionForceAndHessian\n");
	//printf("p: ");
	//CuMatrix::printFloatVec(p, 3);
	//printf("n: ");
	//CuMatrix::printFloatVec(n, 3);
	//printf("c: ");
	//CuMatrix::printFloatVec(c, 3);
	//printf("force: ");
	//CuMatrix::printFloatVec(force, 3);
	//printf("h: ");
	//CuMatrix::printMat3(hessian);


	return;
}


GPU_CPU_INLINE_FUNC  void  GAIA::evaluateNeoHookeanMaterialForceAndHessian(CFloatingTypeGPU miu, CFloatingTypeGPU lmbd, CFloatingTypeGPU A, int32_t tetId, int32_t vOrderInTet,
	const int32_t* tetVIds, CFloatingTypeGPU* vs, CFloatingTypeGPU* vsPrev, CFloatingTypeGPU* DmInvs, FloatingTypeGPU* force, FloatingTypeGPU* h,
	CFloatingTypeGPU dampingVolume, CFloatingTypeGPU dampingShear, CFloatingTypeGPU dt)
{

	CFloatingTypeGPU* v0 = vs + VERTEX_BUFFER_STRIDE * tetVIds[0];
	CFloatingTypeGPU* v1 = vs + VERTEX_BUFFER_STRIDE * tetVIds[1];
	CFloatingTypeGPU* v2 = vs + VERTEX_BUFFER_STRIDE * tetVIds[2];
	CFloatingTypeGPU* v3 = vs + VERTEX_BUFFER_STRIDE * tetVIds[3];

	FloatingTypeGPU Ds[9];
	FloatingTypeGPU F[9];

	CuMatrix::vec3Minus(v1, v0, Ds);
	CuMatrix::vec3Minus(v2, v0, Ds + 3);
	CuMatrix::vec3Minus(v3, v0, Ds + 6);
	CFloatingTypeGPU* DmInv = DmInvs + tetId * 9;
	//printf("Ds:\n");
	//CuMatrix::printMat3(Ds);

	CuMatrix::mat3MatProduct(Ds, DmInv, F);
	//printf("F:\n");
	//CuMatrix::printMat3(F);

	CFloatingTypeGPU a = 1 + miu / lmbd;

	//CFloatingTypeGPU Phi_D = 0.5f * (CuMatrix::mat3FNormSquare(F) - 3.f);
	CFloatingTypeGPU detF = CuMatrix::mat3Determinant(F);
	//CFloatingTypeGPU Phi_H = 0.5f * SQR(detF - a);

	// CFloatingTypeGPU E = A * (miu * Phi_D + lmbd * Phi_H);

	// printf("Phi_D: %f | Phi_H: %f | E: %f\n", Phi_D, detF, Phi_H);

	FloatingTypeGPU* dPhi_D_dF = F;

	// printf("dPhi_D_dF:\n");
	// CuMatrix::printFloatVec(dPhi_D_dF, 9);

	CFloatingTypeGPU F1_1 = F[0];
	CFloatingTypeGPU F2_1 = F[1];
	CFloatingTypeGPU F3_1 = F[2];
	CFloatingTypeGPU F1_2 = F[3 * 1];
	CFloatingTypeGPU F2_2 = F[1 + 3 * 1];
	CFloatingTypeGPU F3_2 = F[2 + 3 * 1];
	CFloatingTypeGPU F1_3 = F[3 * 2];
	CFloatingTypeGPU F2_3 = F[1 + 3 * 2];
	CFloatingTypeGPU F3_3 = F[2 + 3 * 2];

	FloatingTypeGPU ddetF_dF[9] =
	{ F2_2 * F3_3 - F2_3 * F3_2,
	  F1_3 * F3_2 - F1_2 * F3_3,
	  F1_2 * F2_3 - F1_3 * F2_2,
	  F2_3 * F3_1 - F2_1 * F3_3,
	  F1_1 * F3_3 - F1_3 * F3_1,
	  F1_3 * F2_1 - F1_1 * F2_3,
	  F2_1 * F3_2 - F2_2 * F3_1,
	  F1_2 * F3_1 - F1_1 * F3_2,
	  F1_1 * F2_2 - F1_2 * F2_1 };

	// printf("ddetF_dF:\n");
	// CuMatrix::printFloatVec(ddetF_dF, 9);

	CuMatrix::Mat9x9Static<FloatingTypeGPU> d2E_dF_dF;

	CuMatrix::vec9OuterProduct(ddetF_dF, ddetF_dF, d2E_dF_dF);

	CFloatingTypeGPU k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	d2E_dF_dF.multiplyBy(lmbd);

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF.multiplyBy(A);

	// printf("d2E_dF_dF:\n");
	// CuMatrix::printMat(d2E_dF_dF.data, 9, 9);

	//std::cout << "d2E_dF_dF:\n" << d2E_dF_dF << std::endl;

	CuMatrix::vec9Mul(dPhi_D_dF, miu, dPhi_D_dF);
	CuMatrix::vec9Mul(ddetF_dF, lmbd * k, ddetF_dF);

	FloatingTypeGPU dE_dF[9];

	CuMatrix::vec9Add(dPhi_D_dF, ddetF_dF, dE_dF);

	CuMatrix::vec9Mul(dE_dF, A, dE_dF);

	/*printf("dE_dF of tet %d:\n", tetId);
	CuMatrix::printFloatVec(dE_dF, 9);*/

	//std::cout << "dE_dF:\n" << dE_dF << std::endl;

	CFloatingTypeGPU DmInv1_1 = DmInv[0];
	CFloatingTypeGPU DmInv2_1 = DmInv[1];
	CFloatingTypeGPU DmInv3_1 = DmInv[2];
	CFloatingTypeGPU DmInv1_2 = DmInv[3 * 1];
	CFloatingTypeGPU DmInv2_2 = DmInv[1 + 3 * 1];
	CFloatingTypeGPU DmInv3_2 = DmInv[2 + 3 * 1];
	CFloatingTypeGPU DmInv1_3 = DmInv[3 * 2];
	CFloatingTypeGPU DmInv2_3 = DmInv[1 + 3 * 2];
	CFloatingTypeGPU DmInv3_3 = DmInv[2 + 3 * 2];

	CFloatingTypeGPU ms[4][3] = {
		{-DmInv1_1 - DmInv2_1 - DmInv3_1, -DmInv1_2 - DmInv2_2 - DmInv3_2,  -DmInv1_3 - DmInv2_3 - DmInv3_3},
		{DmInv1_1, DmInv1_2, DmInv1_3},
		{DmInv2_1, DmInv2_2, DmInv2_3},
		{DmInv3_1, DmInv3_2, DmInv3_3}
	};

	CFloatingTypeGPU m1 = ms[vOrderInTet][0], m2 = ms[vOrderInTet][1], m3 = ms[vOrderInTet][2];
	assembleVertexVForceAndHessian(dE_dF, &d2E_dF_dF, m1, m2, m3, force, h);

	if (dampingVolume > 0.f || dampingShear > 0.f)
	{
		FloatingTypeGPU dampingH[9];
		CuMatrix::vec9Mul(h, dampingVolume, dampingH);

		FloatingTypeGPU tmp = (m1 * m1 + m2 * m2 + m3 * m3) * miu * A;

		h[0] += tmp;
		h[4] += tmp;
		h[8] += tmp;
		tmp *= dampingShear;
		dampingH[0] += tmp;
		dampingH[4] += tmp;
		dampingH[8] += tmp;
		CuMatrix::vec9Mul(dampingH, 1 / dt, dampingH);

		CFloatingTypeGPU* v = vs + VERTEX_BUFFER_STRIDE * tetVIds[vOrderInTet];
		CFloatingTypeGPU* vPrev = vsPrev + VERTEX_BUFFER_STRIDE * tetVIds[vOrderInTet];
		FloatingTypeGPU displacement[3];
		CuMatrix::vec3Minus(vPrev, v, displacement);

		FloatingTypeGPU dampingForce[3];
		CuMatrix::mat3VecProduct(dampingH, displacement, dampingForce);
		CuMatrix::vec3Add(force, dampingForce, force);

		CuMatrix::vec9Add(h, dampingH, h);
	}

	//printf("dE_dF of tet %d:\n", tetId);
	//CuMatrix::printFloatVec(dE_dF, 9);

	//printf("force of tet %d:\n-------------\n", tetId);
	//CuMatrix::printFloatVec(force, 3); 

}

GPU_CPU_INLINE_FUNC  void GAIA::evaluateNeoHookeanMaterialForceAndDiagonalHessian(FloatingTypeGPU miu, FloatingTypeGPU lmbd, FloatingTypeGPU A, int32_t tetId, int32_t vOrderInTet, const int32_t* tetVIds, FloatingTypeGPU* vs, CFloatingTypeGPU* DmInvs, FloatingTypeGPU* force, FloatingTypeGPU* h)
{
	CFloatingTypeGPU* v0 = vs + 3 * tetVIds[0];
	CFloatingTypeGPU* v1 = vs + 3 * tetVIds[1];
	CFloatingTypeGPU* v2 = vs + 3 * tetVIds[2];
	CFloatingTypeGPU* v3 = vs + 3 * tetVIds[3];

	FloatingTypeGPU Ds[9];
	FloatingTypeGPU F[9];

	CuMatrix::vec3Minus(v1, v0, Ds);
	CuMatrix::vec3Minus(v2, v0, Ds + 3);
	CuMatrix::vec3Minus(v3, v0, Ds + 6);
	CFloatingTypeGPU* DmInv = DmInvs + tetId * 9;
	//printf("Ds:\n");
	//CuMatrix::printMat3(Ds);

	CuMatrix::mat3MatProduct(Ds, DmInv, F);
	//printf("F:\n");
	//CuMatrix::printMat3(F);

	CFloatingTypeGPU a = 1 + miu / lmbd;

	//CFloatingTypeGPU Phi_D = 0.5f * (CuMatrix::mat3FNormSquare(F) - 3.f);
	CFloatingTypeGPU detF = CuMatrix::mat3Determinant(F);
	//CFloatingTypeGPU Phi_H = 0.5f * SQR(detF - a);

	// CFloatingTypeGPU E = A * (miu * Phi_D + lmbd * Phi_H);

	// printf("Phi_D: %f | Phi_H: %f | E: %f\n", Phi_D, detF, Phi_H);

	FloatingTypeGPU* dPhi_D_dF = F;

	// printf("dPhi_D_dF:\n");
	// CuMatrix::printFloatVec(dPhi_D_dF, 9);

	CFloatingTypeGPU F1_1 = F[0];
	CFloatingTypeGPU F2_1 = F[1];
	CFloatingTypeGPU F3_1 = F[2];
	CFloatingTypeGPU F1_2 = F[3 * 1];
	CFloatingTypeGPU F2_2 = F[1 + 3 * 1];
	CFloatingTypeGPU F3_2 = F[2 + 3 * 1];
	CFloatingTypeGPU F1_3 = F[3 * 2];
	CFloatingTypeGPU F2_3 = F[1 + 3 * 2];
	CFloatingTypeGPU F3_3 = F[2 + 3 * 2];

	FloatingTypeGPU ddetF_dF[9] =
	{ F2_2 * F3_3 - F2_3 * F3_2,
	  F1_3 * F3_2 - F1_2 * F3_3,
	  F1_2 * F2_3 - F1_3 * F2_2,
	  F2_3 * F3_1 - F2_1 * F3_3,
	  F1_1 * F3_3 - F1_3 * F3_1,
	  F1_3 * F2_1 - F1_1 * F2_3,
	  F2_1 * F3_2 - F2_2 * F3_1,
	  F1_2 * F3_1 - F1_1 * F3_2,
	  F1_1 * F2_2 - F1_2 * F2_1 };

	// printf("ddetF_dF:\n");
	// CuMatrix::printFloatVec(ddetF_dF, 9);

	CuMatrix::Mat9x9Static<FloatingTypeGPU> d2E_dF_dF;

	CuMatrix::vec9OuterProduct(ddetF_dF, ddetF_dF, d2E_dF_dF);

	CFloatingTypeGPU k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	d2E_dF_dF.multiplyBy(lmbd);

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF.multiplyBy(A);

	// printf("d2E_dF_dF:\n");
	// CuMatrix::printMat(d2E_dF_dF.data, 9, 9);

	//std::cout << "d2E_dF_dF:\n" << d2E_dF_dF << std::endl;

	CuMatrix::vec9Mul(dPhi_D_dF, miu, dPhi_D_dF);
	CuMatrix::vec9Mul(ddetF_dF, lmbd * k, ddetF_dF);

	FloatingTypeGPU dE_dF[9];

	CuMatrix::vec9Add(dPhi_D_dF, ddetF_dF, dE_dF);

	CuMatrix::vec9Mul(dE_dF, A, dE_dF);

	/*printf("dE_dF of tet %d:\n", tetId);
	CuMatrix::printFloatVec(dE_dF, 9);*/

	//std::cout << "dE_dF:\n" << dE_dF << std::endl;

	CFloatingTypeGPU DmInv1_1 = DmInv[0];
	CFloatingTypeGPU DmInv2_1 = DmInv[1];
	CFloatingTypeGPU DmInv3_1 = DmInv[2];
	CFloatingTypeGPU DmInv1_2 = DmInv[3 * 1];
	CFloatingTypeGPU DmInv2_2 = DmInv[1 + 3 * 1];
	CFloatingTypeGPU DmInv3_2 = DmInv[2 + 3 * 1];
	CFloatingTypeGPU DmInv1_3 = DmInv[3 * 2];
	CFloatingTypeGPU DmInv2_3 = DmInv[1 + 3 * 2];
	CFloatingTypeGPU DmInv3_3 = DmInv[2 + 3 * 2];

	FloatingTypeGPU ms[4][3] = {
		{-DmInv1_1 - DmInv2_1 - DmInv3_1, -DmInv1_2 - DmInv2_2 - DmInv3_2,  -DmInv1_3 - DmInv2_3 - DmInv3_3},
		{DmInv1_1, DmInv1_2, DmInv1_3},
		{DmInv2_1, DmInv2_2, DmInv2_3},
		{DmInv3_1, DmInv3_2, DmInv3_3}
	};

	FloatingTypeGPU m1 = ms[vOrderInTet][0];
	FloatingTypeGPU m2 = ms[vOrderInTet][1];
	FloatingTypeGPU m3 = ms[vOrderInTet][2];

	CFloatingTypeGPU A1 = dE_dF[0];
	CFloatingTypeGPU A2 = dE_dF[1];
	CFloatingTypeGPU A3 = dE_dF[2];
	CFloatingTypeGPU A4 = dE_dF[3];
	CFloatingTypeGPU A5 = dE_dF[4];
	CFloatingTypeGPU A6 = dE_dF[5];
	CFloatingTypeGPU A7 = dE_dF[6];
	CFloatingTypeGPU A8 = dE_dF[7];
	CFloatingTypeGPU A9 = dE_dF[8];

	// force is the negative of gradient
	force[0] -= A1 * m1 + A4 * m2 + A7 * m3;
	force[1] -= A2 * m1 + A5 * m2 + A8 * m3;
	force[2] -= A3 * m1 + A6 * m2 + A9 * m3;

	FloatingTypeGPU d2E_dF_dF_X_dF_dXi[9];

	FloatingTypeGPU HL_Row1[9];
	FloatingTypeGPU HL_Row2[9];
	FloatingTypeGPU HL_Row3[9];

	FloatingTypeGPU* HL[3] = { HL_Row1, HL_Row2, HL_Row3 };

	for (int32_t iCol = 0; iCol < 9; iCol++)
	{
		HL_Row1[iCol] = d2E_dF_dF(0, iCol) * m1;
		HL_Row1[iCol] += d2E_dF_dF(3, iCol) * m2;
		HL_Row1[iCol] += d2E_dF_dF(6, iCol) * m3;

		HL_Row2[iCol] = d2E_dF_dF(1, iCol) * m1;
		HL_Row2[iCol] += d2E_dF_dF(4, iCol) * m2;
		HL_Row2[iCol] += d2E_dF_dF(7, iCol) * m3;

		HL_Row3[iCol] = d2E_dF_dF(2, iCol) * m1;
		HL_Row3[iCol] += d2E_dF_dF(5, iCol) * m2;
		HL_Row3[iCol] += d2E_dF_dF(8, iCol) * m3;
	}

	//for (int32_t iRow = 0; iRow < 3; iRow++) {
	//	h[iRow] += HL[iRow][0] * m1 + HL[iRow][3] * m2 + HL[iRow][6] * m3;
	//	h[iRow + 3] += HL[iRow][1] * m1 + HL[iRow][4] * m2 + HL[iRow][7] * m3;
	//	h[iRow + 6] += HL[iRow][2] * m1 + HL[iRow][5] * m2 + HL[iRow][8] * m3;
	//}
	h[0] += HL[0][0] * m1 + HL[0][3] * m2 + HL[0][6] * m3; // H(0,0)
	h[1] += HL[1][1] * m1 + HL[1][4] * m2 + HL[1][7] * m3; // h(1,1)
	h[2] += HL[2][2] * m1 + HL[2][5] * m2 + HL[2][8] * m3; // h(2,2)

}


GPU_CPU_INLINE_FUNC void GAIA::VBDStep_allInOne(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId)
{
	FloatingTypeGPU force[3] = { 0.f, 0.f, 0.f };
	FloatingTypeGPU h[9] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };;

	FloatingTypeGPU* pos = pTetMeshGPU->vertPos + 3 * vertexId;

	//printf("vertex %d's pos: %f %f %f\n", vertexId, pos[0], pos[1], pos[2]);
	//CuMatrix::printFloatVec(pos, 3);

	accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysics->dtSqrReciprocal, force, h);

#ifdef PRINT_DBG_INFO
	printf("after accumulateInertiaForceAndHessian\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif
	accumulateMaterialForceAndHessianForVertex_NeoHookean(pTetMeshGPU, vertexId, force, h, pPhysics->dt);

#ifdef PRINT_DBG_INFO
	printf("after accumulateMaterialForceAndHessian_NeoHookean\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif
	accumulateBoundaryForceAndHessianGPU(pPhysics, pTetMeshGPU, vertexId, force, h);

	if (CuMatrix::vec3NormSquare(force) > CMP_EPSILON2)
	{
		FloatingTypeGPU descentDirection[3];
		FloatingTypeGPU stepSize = pPhysics->stepSize;
		bool solverSuccess = CuMatrix::solve3x3_psd_stable(h, force, descentDirection);

		if (!solverSuccess)
		{
			stepSize = pPhysics->stepSizeGD;
		}

#ifdef PRINT_DBG_INFO
		printf("descentDirection: ");
		CuMatrix::printFloatVec(descentDirection, 3);
#endif
		CuMatrix::vec3MulAddTo(descentDirection, stepSize, pos);
	}

	return;
}

__inline__ __device__ void GAIA::VBDStep_vertexSweep(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU acceleratorOmega)
{
	FloatingTypeGPU force[3] = { 0.f, 0.f, 0.f };
	FloatingTypeGPU h[9] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };;

	FloatingTypeGPU* v = pTetMeshGPU->vertPos + 3 * vertexId;

	//printf("vertex %d's pos: %f %f %f\n", vertexId, pos[0], pos[1], pos[2]);
	//CuMatrix::printFloatVec(pos, 3);

	accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysics->dtSqrReciprocal, force, h);

#ifdef PRINT_DBG_INFO
	printf("after accumulateInertiaForceAndHessian\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif
	accumulateMaterialForceAndHessianForVertex_NeoHookean_preCompuated(pTetMeshGPU, vertexId, force, h);

#ifdef __CUDACC__
	__syncthreads();
#endif // 

	accumulateCollisionForceAndHessian(pPhysics, pTetMeshGPU, vertexId, force, h);
	accumulateBoundaryForceAndHessianGPU(pPhysics, pTetMeshGPU, vertexId, force, h);


#ifdef PRINT_DBG_INFO
	printf("after accumulateMaterialForceAndHessian_NeoHookean\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif

	if (CuMatrix::vec3NormSquare(force) > CMP_EPSILON2)
	{
		FloatingTypeGPU descentDirection[3];
		bool solverSuccess = CuMatrix::solve3x3_psd_stable(h, force, descentDirection);

		FloatingTypeGPU stepSize = solverSuccess ? pPhysics->stepSize : pPhysics->stepSizeGD;

#ifdef PRINT_DBG_INFO
		printf("descentDirection: ");
		CuMatrix::printFloatVec(descentDirection, 3);
#endif
		CuMatrix::vec3MulAddTo(descentDirection, stepSize, v);
	}

	return;
}

__inline__ __device__ void GAIA::VBDStep_vertexSweep_V2(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId_, CFloatingTypeGPU acceleratorOmega)
{

}

GPU_CPU_INLINE_FUNC void GAIA::GDStep_vertexSweep_allInOne(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId)
{
	FloatingTypeGPU force[3] = { 0.f, 0.f, 0.f };
	FloatingTypeGPU h[9] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };;

	CFloatingTypeGPU* pos = pTetMeshGPU->vertPos + 3 * vertexId;
	FloatingTypeGPU* dx = pTetMeshGPU->positionsNew + 3 * vertexId;

	//printf("vertex %d's pos: %f %f %f\n", vertexId, pos[0], pos[1], pos[2]);
	//CuMatrix::printFloatVec(pos, 3);

	accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysics->dtSqrReciprocal, force, h);

	FloatingTypeGPU hDiag[3] = { 0.f, 0.f, 0.f };

	hDiag[0] += h[0];
	hDiag[1] += h[4];
	hDiag[2] += h[8];

	const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
	int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[vertexId * 2];
	int32_t numNeiTets = pTopology->vertexNeighborTets_infos[vertexId * 2 + 1];

	for (int32_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++)
	{
		int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + iNeiTet];
		int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
		int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + iNeiTet];

		evaluateNeoHookeanMaterialForceAndDiagonalHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
			tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->DmInvs, force, hDiag);
	}


#ifdef PRINT_DBG_INFO
	printf("after accumulateInertiaForceAndHessian\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif
	

#ifdef PRINT_DBG_INFO
	printf("after accumulateMaterialForceAndHessian_NeoHookean\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif

	if (CuMatrix::vec3NormSquare(force) > CMP_EPSILON2
		&& hDiag[0] != 0.f
		&& hDiag[1] != 0.f
		&& hDiag[2] != 0.f
		)
	{
		FloatingTypeGPU descentDirection[3] = {
			force[0] / hDiag[0],
			force[1] / hDiag[1],
			force[2] / hDiag[2],
		};

		FloatingTypeGPU stepSize = pPhysics->stepSizeGD;

#ifdef PRINT_DBG_INFO
		printf("descentDirection: ");
		CuMatrix::printFloatVec(descentDirection, 3);
#endif
#ifdef GPU_JACOBI_DX
		CuMatrix::vec3Set(dx, descentDirection);
#endif // GPU_JACOBI_DX

	}
	else
	{
		CuMatrix::vec3Set(dx, 0.f);

	}

	return;
}

GPU_CPU_INLINE_FUNC void GAIA::GDStep_vertexSweep_allInOne_blockJacobi(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId)
{
	FloatingTypeGPU force[3] = { 0.f, 0.f, 0.f };
	FloatingTypeGPU h[9] = { 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f };;

	CFloatingTypeGPU* pos = pTetMeshGPU->vertPos + 3 * vertexId;
	FloatingTypeGPU* dx = pTetMeshGPU->positionsNew + 3 * vertexId;

	//printf("vertex %d's pos: %f %f %f\n", vertexId, pos[0], pos[1], pos[2]);
	//CuMatrix::printFloatVec(pos, 3);

	accumulateInertiaForceAndHessian(pTetMeshGPU, vertexId, pPhysics->dtSqrReciprocal, force, h);


	const TetMeshTopologyGPU* pTopology = pTetMeshGPU->pTopology;
	int32_t neiTetsStart = pTopology->vertexNeighborTets_infos[vertexId * 2];
	int32_t numNeiTets = pTopology->vertexNeighborTets_infos[vertexId * 2 + 1];

	for (int32_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++)
	{
		int32_t neiTetId = pTopology->vertexNeighborTets[neiTetsStart + iNeiTet];
		int32_t* tetVIds = pTopology->tetVIds + 4 * neiTetId;
		int32_t vertexOrderInTet = pTopology->vertexNeighborTets_vertexOrder[neiTetsStart + iNeiTet];

		evaluateNeoHookeanMaterialForceAndHessian(pTetMeshGPU->miu, pTetMeshGPU->lmbd, pTetMeshGPU->tetRestVolume[neiTetId], neiTetId, vertexOrderInTet,
			tetVIds, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->DmInvs, force, h, pTetMeshGPU->dampingVolume, pTetMeshGPU->dampingShear, pPhysics->dt);
	}
	accumulateBoundaryForceAndHessianGPU(pPhysics, pTetMeshGPU, vertexId, force, h);


#ifdef PRINT_DBG_INFO
	printf("after accumulateInertiaForceAndHessian\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif


#ifdef PRINT_DBG_INFO
	printf("after accumulateMaterialForceAndHessian_NeoHookean\n");
	printf("force: ");
	CuMatrix::printFloatVec(force, 3);
	printf("h: ");
	CuMatrix::printMat3(h);
#endif

	if (CuMatrix::vec3NormSquare(force) > CMP_EPSILON2
		)
	{
		FloatingTypeGPU descentDirection[3];
		bool solverSuccess = CuMatrix::solve3x3_psd_stable(h, force, descentDirection);


#ifdef PRINT_DBG_INFO
		printf("descentDirection: ");
		CuMatrix::printFloatVec(descentDirection, 3);
#endif
#ifdef GPU_JACOBI_DX
		if (solverSuccess)
		{
			CuMatrix::vec3Set(dx, descentDirection);
		}
		else
		{
			CuMatrix::vec3Set(dx, 0.f);
		}
#endif // GPU_JACOBI_DX

	}
	else
	{
		CuMatrix::vec3Set(dx, 0.f);

	}

	return;
}



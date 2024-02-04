#pragma once

#include "PBDTetMeshNeoHookeanCompute.h"
#include "../Parallelization//GPUParallelization.h"
#include "../Parallelization//CPUParallelization.h"
#include "PBDTetMeshGeneralCompute.h"

//#define AGGREGATED_SOLVE

using namespace GAIA;
  
__device__ void printGPUMesh_NeoHookean(GAIA::TetMeshFEMGPU_NeoHookean* pTetMesh) {
	printf("nVerts: %d\n", pTetMesh->nVerts);
	printf("nTets: %d\n", pTetMesh->nTets);
	printf("tetRestVolume: %p\n", pTetMesh->tetRestVolume);
	printf("tetInvRestVolume: %p\n", pTetMesh->tetInvRestVolume);
	printf("vertexInvMass: %p\n", pTetMesh->vertexInvMass);
	printf("DSInvs: %p\n", pTetMesh->DmInvs);
	printf("tetVIds: %p\n", pTetMesh->tetVIds);
	printf("vertPos: %p\n", pTetMesh->vertPos);
	printf("vertPrevPos: %p\n", pTetMesh->vertPrevPos);
	printf("velocity: %p\n", pTetMesh->velocity);

	printf("devCompliance: %f\n", pTetMesh->devCompliance);
	printf("devDamping: %f\n", pTetMesh->devDamping);
	printf("volCompliance: %f\n", pTetMesh->volCompliance);

	printf("devLambdas: %p\n", pTetMesh->devLambdas);
	printf("volLambdas: %p\n", pTetMesh->volLambdas);
}

void GAIA::solveMaterial(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	int32_t numThreads, TetMeshFEMGPU_NeoHookean* pTetMeshGPU, cudaStream_t stream)
{
	//parallel_for_tets KERNEL_ARGS4((tetParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream) 
	//	(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU, solveMaterialForOneTet);

	// v1
	if (tetParallelizationGroupSize > numThreads)
	{
		solveMaterialGPU KERNEL_ARGS4((tetParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU);
	}
	else
	{
		solveMaterialGPU KERNEL_ARGS4(1, tetParallelizationGroupSize, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU);
	}

	// v2
	//if (tetParallelizationGroupSize > numThreads )
	//{
	//	solveMaterialGPU KERNEL_ARGS4((tetParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)
	//		(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->tetVIds, pTetMeshGPU->vertexInvMass,
	//			pTetMeshGPU->tetInvRestVolume, pTetMeshGPU->DmInvs, pTetMeshGPU->devLambdas, pTetMeshGPU->dt,
	//			pTetMeshGPU->devCompliance, pTetMeshGPU->devDamping);
	//}
	//else
	//{
	//	solveMaterialGPU KERNEL_ARGS4(1, tetParallelizationGroupSize, 0, stream)
	//		(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->tetVIds, pTetMeshGPU->vertexInvMass,
	//			pTetMeshGPU->tetInvRestVolume, pTetMeshGPU->DmInvs, pTetMeshGPU->devLambdas, pTetMeshGPU->dt,
	//			pTetMeshGPU->devCompliance, pTetMeshGPU->devDamping);
	//}

}

void GAIA::solveMaterial_oneThreadForLoop(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, int32_t numThreads, 
	TetMeshFEMGPU_NeoHookean* pTetMeshGPU, cudaStream_t stream)
{
	//printf("Size of TetMeshFEMGPU_NeoHookean on CPU (.cu): %d\n", sizeof(TetMeshFEMGPU_NeoHookean));
	// serial: for debug only: 
	solveMaterialGPU_oneThreadForLoop KERNEL_ARGS4(1,1,0, stream) (*pTetMeshGPU);
}

__global__ void GAIA::solveMaterialGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, int32_t* tetVIdsAll, FloatingTypeGPU* vertexInvMass,
	FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU* DmInvs, FloatingTypeGPU* devLambdas, FloatingTypeGPU dt, 
	FloatingTypeGPU devCompliance,	FloatingTypeGPU devDamping)
{
	int iTet = blockIdx.x * blockDim.x + threadIdx.x;

	if (iTet < tetParallelizationGroupSize)
	{
		int32_t tetId = pTetParallelizationGroup[iTet];
		int32_t* tetVIds = tetVIdsAll + tetId * 4;
		int32_t id0 = tetVIds[0];
		int32_t id1 = tetVIds[1];
		int32_t id2 = tetVIds[2];
		int32_t id3 = tetVIds[3];

#ifdef EVAL_TET_ACCESS_COUNT
		pTetMeshGPU->tetAccessCount[tetId] += 1;
		pTetMeshGPU->vertAccessCount[id0] += 1;
		pTetMeshGPU->vertAccessCount[id1] += 1;
		pTetMeshGPU->vertAccessCount[id2] += 1;
		pTetMeshGPU->vertAccessCount[id3] += 1;
#endif // EVAL_TET_ACCESS_COUNT

		FloatingTypeGPU* v0 = vertPos + id0 * 3;
		FloatingTypeGPU* v1 = vertPos + id1 * 3;
		FloatingTypeGPU* v2 = vertPos + id2 * 3;
		FloatingTypeGPU* v3 = vertPos + id3 * 3;

		FloatingTypeGPU gradDev[12];
		FloatingTypeGPU gradVol[12];

		FloatingTypeGPU CDev = solveDevConstraint(tetId, v0, v1, v2, v3, gradDev, DmInvs);
		applyToElem(vertPos, vertPrevPos, vertexInvMass, tetInvRestVolume, tetId, tetVIds, CDev, devCompliance, dt,
			devDamping, devLambdas + tetId, gradDev);

		FloatingTypeGPU CVol = solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, DmInvs);
		applyToInfiniteStiffness(vertPos, vertexInvMass, tetId, tetVIds, CVol, gradVol);
	}
	return ;
}

__global__ void GAIA::solveMaterialGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	int iTet = blockIdx.x * blockDim.x + threadIdx.x;
	if (iTet < tetParallelizationGroupSize)
	{
		int32_t tetId = pTetParallelizationGroup[iTet];
		solveMaterialForOneTet(tetId, pTetMeshGPU);
	}
	return;
}

__global__ void GAIA::solveNeoHookeanMaterialAggregatedGPU_kernel(TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	int iThread = blockIdx.x * blockDim.x + threadIdx.x;
	for (size_t iParalGroup = 0; iParalGroup < pTetMeshGPU->numTetsColoredCatergories; iParalGroup++)
	{
		if (iThread >= pTetMeshGPU->tetsColoringEachCategorySize[iParalGroup]) {
			return;
		}

		int32_t tetId = pTetMeshGPU->tetsColoringCategoriesPointers[iParalGroup][iThread];
		solveMaterialForOneTet(tetId, pTetMeshGPU);
	}
}

void GAIA::solveNeoHookeanMaterialAggregatedGPU(int32_t numThreads, cudaStream_t stream, int32_t maxParallelGroupSize, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	solveNeoHookeanMaterialAggregatedGPU_kernel KERNEL_ARGS4((maxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)	(pTetMeshGPU);
}

__global__ void GAIA::solveMaterialGPU_oneThreadForLoop(TetMeshFEMGPU_NeoHookean tetMeshGPU)
{
	//printf("Size of TetMeshFEMGPU_NeoHookean on GPU: %d\n", sizeof(TetMeshFEMGPU_NeoHookean));
	//printGPUMesh_NeoHookean(&tetMeshGPU);
	//return;

	for (int iTet = 0; iTet < tetMeshGPU.nTets; iTet++)
	{
		solveMaterialForOneTet(iTet, &tetMeshGPU);
	}
}


void GAIA::solveMaterialParallelOnCPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	cpu_parallel_for(0, tetParallelizationGroupSize, [&](int iTet)
		{
			int32_t tetId = pTetParallelizationGroup[iTet];
			const int32_t* tetVIds = pTetMeshGPU->tetVIds + tetId * 4;

			int32_t id0 = tetVIds[0];
			int32_t id1 = tetVIds[1];
			int32_t id2 = tetVIds[2];
			int32_t id3 = tetVIds[3];

			FloatingTypeGPU* v0 = pTetMeshGPU->vertPos + id0 * 3;
			FloatingTypeGPU* v1 = pTetMeshGPU->vertPos + id1 * 3;
			FloatingTypeGPU* v2 = pTetMeshGPU->vertPos + id2 * 3;
			FloatingTypeGPU* v3 = pTetMeshGPU->vertPos + id3 * 3;

			 //printf("TId: %d\nv0: %f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
			 //	tetId,
			 //	v0[0], v0[1], v0[1],
			 //	v1[0], v1[1], v1[1],
			 //	v2[0], v2[1], v2[1],
			 //	v3[0], v3[1], v3[1]);

			FloatingTypeGPU gradDev[12];
			FloatingTypeGPU gradVol[12];

			TetMeshFEMGPU_NeoHookean* pTetMeshGPU_NeoHookean = (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU;
			FloatingTypeGPU CDev = solveDevConstraint(tetId, v0, v1, v2, v3, gradDev, pTetMeshGPU_NeoHookean->DmInvs);
			applyToElem(pTetMeshGPU, tetId, tetVIds, CDev, pTetMeshGPU_NeoHookean->devCompliance, pTetMeshGPU->dt,
				pTetMeshGPU_NeoHookean->devDamping, pTetMeshGPU_NeoHookean->devLambdas + tetId, gradDev);

			FloatingTypeGPU CVol = solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, pTetMeshGPU->DmInvs);
			applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, tetId, tetVIds, CVol, gradVol);
		});
}

void GAIA::test_MaterialSolve_oneTet_call(int32_t iTet, TetMeshFEMGPU_NeoHookean* pTetMeshGPU, FloatingTypeGPU* gradDev, FloatingTypeGPU* CDev)
{
	test_solveDev_oneTet KERNEL_ARGS2(1, 1) (iTet, pTetMeshGPU->vertPos, pTetMeshGPU->vertPrevPos, pTetMeshGPU->tetVIds, pTetMeshGPU->vertexInvMass,
		pTetMeshGPU->tetInvRestVolume, pTetMeshGPU->DmInvs, pTetMeshGPU->devLambdas, pTetMeshGPU->dt,
		pTetMeshGPU->devCompliance, pTetMeshGPU->devDamping, gradDev, CDev);
}

__global__ void GAIA::test_solveDev_oneTet(int32_t tetId, FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos,
	const int32_t* tetVIdsAll, FloatingTypeGPU* vertexInvMass, FloatingTypeGPU* tetInvRestVolume,
	FloatingTypeGPU* DmInvs, FloatingTypeGPU* devLambdas, FloatingTypeGPU dt, FloatingTypeGPU devCompliance,
	FloatingTypeGPU devDamping, FloatingTypeGPU* gradDev, FloatingTypeGPU* CDev)
{
	const int32_t* tetVIds = tetVIdsAll + tetId * 4;
	int32_t id0 = tetVIds[0];
	int32_t id1 = tetVIds[1];
	int32_t id2 = tetVIds[2];
	int32_t id3 = tetVIds[3];

	FloatingTypeGPU* v0 = vertPos + id0 * 3;
	FloatingTypeGPU* v1 = vertPos + id1 * 3;
	FloatingTypeGPU* v2 = vertPos + id2 * 3;
	FloatingTypeGPU* v3 = vertPos + id3 * 3;

	//printf("TId: %d\nv0: %f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
	//tetId,
	//v0[0], v0[1], v0[1],
	//v1[0], v1[1], v1[1],
	//v2[0], v2[1], v2[1],
	//v3[0], v3[1], v3[1]);

	*CDev = solveDevConstraint(tetId, v0, v1, v2, v3, gradDev, DmInvs);
	//printf("gradDev (before apply to elem): ");
	//for (size_t i = 0; i < 12; i++)
	//{
	//	printf(" %e", gradDev[i]);
	//}
	//printf("\ntest_MaterialSolve_oneTet on gpu done!\n");

	applyToElem(vertPos, vertPrevPos, vertexInvMass, tetInvRestVolume, tetId, tetVIds, *CDev, devCompliance, dt,
		devDamping, devLambdas + tetId, gradDev);

}

__host__ __device__ void GAIA::solveMaterialForOneTet(int32_t tetId, FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, int32_t* tetVIds, FloatingTypeGPU* vertexInvMass, FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU* DmInvs)
{
	return ;
}

__forceinline __host__ __device__ void GAIA::solveMaterialForOneTet(int32_t tetId, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	//printf("Solving material for: %d\n", tetId);

	const int32_t* tetVIds = pTetMeshGPU->tetVIds + tetId * 4;
	int32_t id0 = tetVIds[0];
	int32_t id1 = tetVIds[1];
	int32_t id2 = tetVIds[2];
	int32_t id3 = tetVIds[3];

#ifdef EVAL_TET_ACCESS_COUNT
	pTetMeshGPU->tetAccessCount[tetId] += 1;
	pTetMeshGPU->vertAccessCount[id0] += 1;
	pTetMeshGPU->vertAccessCount[id1] += 1;
	pTetMeshGPU->vertAccessCount[id2] += 1;
	pTetMeshGPU->vertAccessCount[id3] += 1;
#endif // EVAL_TET_ACCESS_COUNT

	//if (tetId > pTetMeshGPU->nTets)
	//{
	//	printf("!!! tetId: %d > pTetMeshGPU->nTets: %d\n", tetId, pTetMeshGPU->nTets);
	//}

	//for (int i = 0; i < 4; i++)
	//{
	//	if (tetVIds[i] > pTetMeshGPU->nVerts)
	//	{
	//		printf("!!! vid: %d > pTetMeshGPU->nVerts: %d\n", tetVIds[0], pTetMeshGPU->nVerts);
	//	}
	//}

	FloatingTypeGPU* v0 = pTetMeshGPU->vertPos + id0 * 3;
	FloatingTypeGPU* v1 = pTetMeshGPU->vertPos + id1 * 3;
	FloatingTypeGPU* v2 = pTetMeshGPU->vertPos + id2 * 3;
	FloatingTypeGPU* v3 = pTetMeshGPU->vertPos + id3 * 3;

	// printf("TId: %d\nv0: %f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
	// 	tetId,
	// 	v0[0], v0[1], v0[1],
	// 	v1[0], v1[1], v1[1],
	// 	v2[0], v2[1], v2[1],
	// 	v3[0], v3[1], v3[1]);

	FloatingTypeGPU gradDev[12];
	FloatingTypeGPU gradVol[12];

#ifdef AGGREGATED_SOLVE
	FloatingTypeGPU CDev = solveDevConstraint(tetId, v0, v1, v2, v3, gradDev, pTetMeshGPU);

	FloatingTypeGPU CVol = solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU);
	
	applyToElem(pTetMeshGPU, tetId, tetVIds, CDev, pTetMeshGPU->devCompliance, pTetMeshGPU->dt,
		pTetMeshGPU->devDamping, pTetMeshGPU->devLambdas + tetId, gradDev);

	if (pTetMeshGPU->volCompliance = 0.f)
	{
		applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, tetId, tetVIds, CVol, gradVol);
	}
	else {
		applyToElem(pTetMeshGPU, tetId, tetVIds, CVol, pTetMeshGPU->volCompliance, pTetMeshGPU->dt,
			pTetMeshGPU->volCompliance, pTetMeshGPU->volLambdas + tetId, gradVol);
	}
#else
	FloatingTypeGPU CDev = solveDevConstraint(tetId, v0, v1, v2, v3, gradDev, pTetMeshGPU);
	applyToElem(pTetMeshGPU, tetId, tetVIds, CDev, pTetMeshGPU->devCompliance, pTetMeshGPU->dt,
		pTetMeshGPU->devDamping, pTetMeshGPU->devLambdas + tetId, gradDev);

	FloatingTypeGPU CVol = solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU);

	if (pTetMeshGPU->volCompliance = 0.f)
	{
		applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, tetId, tetVIds, CVol, gradVol);
	}
	else {
		applyToElem(pTetMeshGPU, tetId, tetVIds, CVol, pTetMeshGPU->volCompliance, pTetMeshGPU->dt,
			pTetMeshGPU->volCompliance, pTetMeshGPU->volLambdas + tetId, gradVol);
	}
#endif // AGGREGATED_SOLVE



	//#ifdef EVAL_TET_ACCESS_COUNT
	//		pTetMeshGPU->vertAccessCount[id0] -= 1;
	//		pTetMeshGPU->vertAccessCount[id1] -= 1;
	//		pTetMeshGPU->vertAccessCount[id2] -= 1;
	//		pTetMeshGPU->vertAccessCount[id3] -= 1;
	//#endif // EVAL_TET_ACCESS_COUNT
//	return ;
}

__host__ __device__ float GAIA::solveDevConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
	FloatingTypeGPU* gradDev, FloatingTypeGPU* DmInvs)
{
	//printf("solveDevConstraint\nv0:\n%f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
	//	v0[0], v0[1], v0[1],
	//	v1[0], v1[1], v1[1],
	//	v2[0], v2[1], v2[1],
	//	v3[0], v3[1], v3[1]);

	FloatingTypeGPU C = 0.0f;
	FloatingTypeGPU P[9];
	FloatingTypeGPU F[9];

	CuMatrix::vec3Minus(v1, v0, P);
	CuMatrix::vec3Minus(v2, v0, P+3);
	CuMatrix::vec3Minus(v3, v0, P+6);
	FloatingTypeGPU* invDS = DmInvs + tetId * 9;
	CuMatrix::mat3MatProduct(P, invDS, F);

	FloatingTypeGPU F0Norm2 = CuMatrix::vec3NormSquare(F);
	FloatingTypeGPU F1Norm2 = CuMatrix::vec3NormSquare(F+3);
	FloatingTypeGPU F2Norm2 = CuMatrix::vec3NormSquare(F+6);

	FloatingTypeGPU r_s = sqrt(F0Norm2 + F1Norm2 + F2Norm2);

	FloatingTypeGPU r_s_inv = r_s != 0.f ? 1.0f / r_s : 0.f;

	FloatingTypeGPU temp[3];

	FloatingTypeGPU* g1 = gradDev + 3;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 0), g1);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 1), temp);
	CuMatrix::vec3Add(temp, g1, g1);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 2), temp);
	CuMatrix::vec3Add(temp, g1, g1);

	FloatingTypeGPU* g2 = gradDev + 6;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 0), g2);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 1), temp);
	CuMatrix::vec3Add(temp, g2, g2);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 2), temp);
	CuMatrix::vec3Add(temp, g2, g2);

	// Wrong at this code block
	FloatingTypeGPU* g3 = gradDev + 9;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 0), g3);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 1), temp);
	CuMatrix::vec3Add(temp, g3, g3);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 2), temp);
	CuMatrix::vec3Add(temp, g3, g3);
	
	gradDev[0] = 0.f;
	gradDev[1] = 0.f;
	gradDev[2] = 0.f;

	CuMatrix::vec3Minus(gradDev, g1, gradDev);
	CuMatrix::vec3Minus(gradDev, g2, gradDev);
	CuMatrix::vec3Minus(gradDev, g3, gradDev);

	C = r_s;

	return C;

	//if (curPhysics.restStableDevProjection) {
	//	// the original dev energy is non-zero
	//	C = r_s - 1.732050807568877;


}

__host__ __device__ float GAIA::solveDevConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, 
	FloatingTypeGPU* v3, FloatingTypeGPU* gradDev, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	FloatingTypeGPU C = 0.0f;
	FloatingTypeGPU P[9];
	FloatingTypeGPU F[9];

	CuMatrix::vec3Minus(v1, v0, P);
	CuMatrix::vec3Minus(v2, v0, P + 3);
	CuMatrix::vec3Minus(v3, v0, P + 6);
	FloatingTypeGPU* invDS = pTetMeshGPU->DmInvs + tetId * 9;
	CuMatrix::mat3MatProduct(P, invDS, F);

	FloatingTypeGPU F0Norm2 = CuMatrix::vec3NormSquare(F);
	FloatingTypeGPU F1Norm2 = CuMatrix::vec3NormSquare(F + 3);
	FloatingTypeGPU F2Norm2 = CuMatrix::vec3NormSquare(F + 6);

	FloatingTypeGPU r_s = sqrtf(F0Norm2 + F1Norm2 + F2Norm2);

	FloatingTypeGPU r_s_inv = r_s != 0.f ? 1.0f / r_s : 0.f;

	FloatingTypeGPU temp[3];

	FloatingTypeGPU* g1 = gradDev + 3;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 0), g1);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 1), temp);
	CuMatrix::vec3Add(temp, g1, g1);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 0, 2), temp);
	CuMatrix::vec3Add(temp, g1, g1);

	FloatingTypeGPU* g2 = gradDev + 6;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 0), g2);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 1), temp);
	CuMatrix::vec3Add(temp, g2, g2);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 1, 2), temp);
	CuMatrix::vec3Add(temp, g2, g2);

	// Wrong at this code block
	FloatingTypeGPU* g3 = gradDev + 9;
	CuMatrix::vec3Mul(F, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 0), g3);
	CuMatrix::vec3Mul(F + 3, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 1), temp);
	CuMatrix::vec3Add(temp, g3, g3);
	CuMatrix::vec3Mul(F + 6, r_s_inv * CuMatrix::mat3IJ(invDS, 2, 2), temp);
	CuMatrix::vec3Add(temp, g3, g3);

	gradDev[0] = 0.f;
	gradDev[1] = 0.f;
	gradDev[2] = 0.f;

	CuMatrix::vec3Minus(gradDev, g1, gradDev);
	CuMatrix::vec3Minus(gradDev, g2, gradDev);
	CuMatrix::vec3Minus(gradDev, g3, gradDev);

	C = r_s;

	return C;
}


__host__ __device__ FloatingTypeGPU GAIA::solveVolConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2,
	FloatingTypeGPU* v3, FloatingTypeGPU* gradVol, TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	FloatingTypeGPU C = 0.0f;
	FloatingTypeGPU P[9];
	FloatingTypeGPU F[9];

	//printf("v0:\n%f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
	//	v0[0], v0[1], v0[1],
	//	v1[0], v1[1], v1[1],
	//	v2[0], v2[1], v2[1],
	//	v3[0], v3[1], v3[1]);


	CuMatrix::vec3Minus(v1, v0, P);
	CuMatrix::vec3Minus(v2, v0, P + 3);
	CuMatrix::vec3Minus(v3, v0, P + 6);
	FloatingTypeGPU* invDS = pTetMeshGPU->DmInvs + tetId * 9;

	CuMatrix::mat3MatProduct(P, invDS, F);

	FloatingTypeGPU* F0 = F;
	FloatingTypeGPU* F1 = F + 3;
	FloatingTypeGPU* F2 = F + 6;

	FloatingTypeGPU dF0[3]; // = ((F1) ^ (F2));
	CuMatrix::vec3CrossProduct(F1, F2, dF0);
	FloatingTypeGPU dF1[3]; // = ((F2) ^ (F0));
	CuMatrix::vec3CrossProduct(F2, F0, dF1);
	FloatingTypeGPU dF2[3]; // = ((F0) ^ (F1));
	CuMatrix::vec3CrossProduct(F0, F1, dF2);

	FloatingTypeGPU temp[3];


	FloatingTypeGPU* g1 = gradVol + 3;
	CuMatrix::vec3Mul(dF0, CuMatrix::mat3IJ(invDS, 0, 0), g1);
	CuMatrix::vec3Mul(dF1, CuMatrix::mat3IJ(invDS, 0, 1), temp);
	CuMatrix::vec3Add(temp, g1, g1);
	CuMatrix::vec3Mul(dF2, CuMatrix::mat3IJ(invDS, 0, 2), temp);
	CuMatrix::vec3Add(temp, g1, g1);

	FloatingTypeGPU* g2 = gradVol + 6;
	CuMatrix::vec3Mul(dF0, CuMatrix::mat3IJ(invDS, 1, 0), g2);
	CuMatrix::vec3Mul(dF1, CuMatrix::mat3IJ(invDS, 1, 1), temp);
	CuMatrix::vec3Add(temp, g2, g2);
	CuMatrix::vec3Mul(dF2, CuMatrix::mat3IJ(invDS, 1, 2), temp);
	CuMatrix::vec3Add(temp, g2, g2);

	FloatingTypeGPU* g3 = gradVol + 9;
	CuMatrix::vec3Mul(dF0, CuMatrix::mat3IJ(invDS, 2, 0), g3);
	CuMatrix::vec3Mul(dF1, CuMatrix::mat3IJ(invDS, 2, 1), temp);
	CuMatrix::vec3Add(temp, g3, g3);
	CuMatrix::vec3Mul(dF2, CuMatrix::mat3IJ(invDS, 2, 2), temp);
	CuMatrix::vec3Add(temp, g3, g3);

	// g0
	gradVol[0] = 0.f;
	gradVol[1] = 0.f;
	gradVol[2] = 0.f;

	CuMatrix::vec3Minus(gradVol, g1, gradVol);
	CuMatrix::vec3Minus(gradVol, g2, gradVol);
	CuMatrix::vec3Minus(gradVol, g3, gradVol);

	FloatingTypeGPU vol = CuMatrix::mat3Determinant(F);

	//printf("F:\n%f %f %f\n%f %f %f\n%f %f %f\ninvDS:\n%f %f %f\n%f %f %f\n%f %f %f\nP:\n%f %f %f\n%f %f %f\n%f %f %f\n",
	//	F[0], F[3], F[6],
	//	F[1], F[4], F[7],
	//	F[2], F[5], F[8],
	//	invDS[0], invDS[3], invDS[6],
	//	invDS[1], invDS[4], invDS[7],
	//	invDS[2], invDS[5], invDS[8],
	//	P[0], P[3], P[6],
	//	P[1], P[4], P[7],
	//	P[2], P[5], P[8]
	//	);

/*	printf("g0:%f %f %f\ng1: %f %f %f\ng2: %f %f %f\ng3: %f %f %f\n",
		gradVol[0], gradVol[1], gradVol[2],
		g1[0], g1[1], g1[2],
		g2[0], g2[1], g2[2],
		g3[0], g3[1], g3[2]
	);	*/


	//C = vol - 1.0 - (mu / lambda);
	C = vol - 1.0f;

	//printf("Vol: %f | Vol loss: %f\n", vol, C);

	return C;
}

__host__ __device__ void GAIA::applyToElem(FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, FloatingTypeGPU* vertexInvMass,
	FloatingTypeGPU* tetInvRestVolume, int32_t tetId, const int32_t* tetVIds,
	FloatingTypeGPU C, FloatingTypeGPU compliance, FloatingTypeGPU dt, FloatingTypeGPU damping, FloatingTypeGPU* constraintLambda, FloatingTypeGPU* gradients)
{
	//printf("damping: %f, constraintLambda: %p, gradients: %p\n", damping, constraintLambda, gradients);

	if (C == 0.0f) {
		return;
	}

	FloatingTypeGPU* gs[4] = {
		gradients,
		gradients + 3,
		gradients + 6,
		gradients + 9
	};

	FloatingTypeGPU w = 0.0f;
	FloatingTypeGPU damp = 0.0f;
	
	bool apply_damp = damping != 0.f;
	
	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		w += CuMatrix::vec3NormSquare(gs[i]) * vertexInvMass[vId];

		if (apply_damp)
		{
			FloatingTypeGPU vel[3];
			CuMatrix::vec3Minus(vertPos + 3 * vId, vertPrevPos + 3 * vId, vel);

			damp += CuMatrix::vec3DotProduct(gs[i], vel);
			// CPoint vel = *(mesh.verts[id]) - mesh.prevPos[id];
			// damp += gradients[i] * vel;
		}
	}


	if (w == 0.0f) {
		return;
	}
	if (apply_damp) {
		C += damp * (damping) / dt;
		w *= 1 + damping / dt;
	}

	FloatingTypeGPU correctedCompliance = tetInvRestVolume[tetId] * compliance / (dt * dt);

	FloatingTypeGPU dlambda = 0.f;
	if (w + correctedCompliance != 0.f) {
		dlambda = (-C - correctedCompliance * (*constraintLambda)) / (w + correctedCompliance);
	}

	*constraintLambda += dlambda;

	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		FloatingTypeGPU* v = vertPos + vId * 3;

		//printf("dlambda: %e,  vertexInvMass[vId]: %e\n", dlambda, vertexInvMass[vId]);

		CuMatrix::vec3Mul(gs[i], dlambda * vertexInvMass[vId], gs[i]);

		//printf("gs[%d]: %e, %e, %e\n", i, gs[i][0], gs[i][1], gs[i][2]);

		CuMatrix::vec3Add(v, gs[i], v);

		//*(mesh.verts[id]) = (*(mesh.verts[id])) + ((dlambda) * (invMass[id]) * (gradients[i]));

	}
}

__host__ __device__ void GAIA::applyToElem(PBDTetMeshFEMGPU* pTetMeshGPU, int32_t tetId, const int32_t* tetVIds, FloatingTypeGPU C,
	FloatingTypeGPU compliance, FloatingTypeGPU dt, FloatingTypeGPU damping, FloatingTypeGPU* constraintLambda, FloatingTypeGPU* gradients)
{
	//printf("damping: %f, constraintLambda: %p, gradients: %p\n", damping, constraintLambda, gradients);

	if (C == 0.0f) {
		return;
	}

	FloatingTypeGPU* gs[4] = {
		gradients,
		gradients + 3,
		gradients + 6,
		gradients + 9
	};

	FloatingTypeGPU w = 0.0f;
	FloatingTypeGPU damp = 0.0f;
	bool apply_damp = damping != 0.f;

	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		w += CuMatrix::vec3NormSquare(gs[i]) * pTetMeshGPU->vertexInvMass[vId];

		if (apply_damp)
		{
			FloatingTypeGPU vel[3];
			CuMatrix::vec3Minus(pTetMeshGPU->vertPos + 3 * vId, pTetMeshGPU->vertPrevPos + 3 * vId, vel);

			damp += CuMatrix::vec3DotProduct(gs[i], vel);
			// CPoint vel = *(mesh.verts[id]) - mesh.prevPos[id];
			// damp += gradients[i] * vel;
		}
	}


	if (w == 0.0f) {
		return;
	}
	if (apply_damp)
	{
		C += damp * (damping) / dt;
		w *= 1 + damping / dt;
	}

	FloatingTypeGPU correctedCompliance = pTetMeshGPU->tetInvRestVolume[tetId] * compliance / (dt * dt);

	FloatingTypeGPU dlambda = 0.f;
	if (w + correctedCompliance != 0.f) {
		dlambda = (-C - correctedCompliance * (*constraintLambda)) / (w + correctedCompliance);
	}

	*constraintLambda += dlambda;

	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		FloatingTypeGPU* v = pTetMeshGPU->vertPos + vId * 3;

		CuMatrix::vec3Mul(gs[i], dlambda * pTetMeshGPU->vertexInvMass[vId], gs[i]);

		CuMatrix::vec3Add(v, gs[i], v);

		//*(mesh.verts[id]) = (*(mesh.verts[id])) + ((dlambda) * (invMass[id]) * (gradients[i]));

	}
}

void GAIA::solveMaterialSequantialOnCPU(TetMeshFEMGPU_NeoHookean* pTetMeshGPU)
{
	for (int32_t iTet = 0; iTet < pTetMeshGPU->nTets; iTet++)
	{
		const int32_t* tetVIds = pTetMeshGPU->tetVIds + iTet * 4;
		int32_t id0 = tetVIds[0];
		int32_t id1 = tetVIds[1];
		int32_t id2 = tetVIds[2];
		int32_t id3 = tetVIds[3];

		FloatingTypeGPU* v0 = pTetMeshGPU->vertPos + id0 * 3;
		FloatingTypeGPU* v1 = pTetMeshGPU->vertPos + id1 * 3;
		FloatingTypeGPU* v2 = pTetMeshGPU->vertPos + id2 * 3;
		FloatingTypeGPU* v3 = pTetMeshGPU->vertPos + id3 * 3;



		// printf("TId: %d\nv0: %f %f %f\nv1: %f %f %f\nv2: %f %f %f\nv3: %f %f %f\n",
		// 	tetId,
		// 	v0[0], v0[1], v0[1],
		// 	v1[0], v1[1], v1[1],
		// 	v2[0], v2[1], v2[1],
		// 	v3[0], v3[1], v3[1]);

		FloatingTypeGPU gradDev[12];
		FloatingTypeGPU gradVol[12];

		TetMeshFEMGPU_NeoHookean* pTetMeshGPU_NeoHookean = (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU;
		FloatingTypeGPU CDev = solveDevConstraint(iTet, v0, v1, v2, v3, gradDev, pTetMeshGPU_NeoHookean->DmInvs);
		applyToElem(pTetMeshGPU, iTet, tetVIds, CDev, pTetMeshGPU_NeoHookean->devCompliance, pTetMeshGPU->dt,
			pTetMeshGPU_NeoHookean->devDamping, pTetMeshGPU_NeoHookean->devLambdas + iTet, gradDev);

		FloatingTypeGPU CVol = solveVolConstraint(iTet, v0, v1, v2, v3, gradVol, pTetMeshGPU->DmInvs);
		applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, iTet, tetVIds, CVol, gradVol);
	}
}


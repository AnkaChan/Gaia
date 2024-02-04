#pragma once

#include <CuMatrix/Buffers/ManagedBuffer.h>
#include "PBDTetMeshFEMGPU.h"
#include <CuMatrix/MatrixOps/CuMatrix.h>


#define COMPILE_TEST_FUNCTIONS

namespace GAIA {
	struct TetMeshFEMGPU_NeoHookean : public PBDTetMeshFEMGPU {
		typedef std::shared_ptr<TetMeshFEMGPU_NeoHookean> SharedPtr;
		typedef TetMeshFEMGPU_NeoHookean* Ptr;

		FloatingTypeGPU devCompliance;
		FloatingTypeGPU devDamping;
		FloatingTypeGPU volCompliance;

		FloatingTypeGPU* devLambdas;
		FloatingTypeGPU* volLambdas;

	};

	__host__ __device__ void solveMaterialForOneTet(int32_t tetId, FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos,
		int32_t* tetVIds, FloatingTypeGPU* vertexInvMass,
		FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU* DmInvs);

	__host__ __device__ void solveMaterialForOneTet(int32_t iTet, TetMeshFEMGPU_NeoHookean* pTetMeshGPU);

	//__host__ __device__ float solveDevConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
	//	FloatingTypeGPU* gradDev, TetMeshFEMGPU_NeoHookean* pTetMeshGPU);

	__host__ __device__ float solveDevConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
		FloatingTypeGPU* gradDev, FloatingTypeGPU* DmInvs);

	__host__ __device__ float solveDevConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
		FloatingTypeGPU* gradDev, TetMeshFEMGPU_NeoHookean* pTetMeshGPU);



	__host__  __device__ FloatingTypeGPU solveVolConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
		FloatingTypeGPU* gradVol, TetMeshFEMGPU_NeoHookean* pTetMeshGPU);

	//__inline__ __host__  __device__ void applyToElem(TetMeshFEMGPU* pTetMeshGPU, int32_t tetId, int32_t* tetVIds,
	//	//FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
	//	FloatingTypeGPU C, FloatingTypeGPU compliance, FloatingTypeGPU dt, FloatingTypeGPU damping,
	//	FloatingTypeGPU* constraintLambda, FloatingTypeGPU* gradients);


	__host__  __device__ void applyToElem(FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, FloatingTypeGPU* vertexInvMass,
		FloatingTypeGPU* tetInvRestVolume, int32_t tetId, const int32_t* tetVIds,
		FloatingTypeGPU C, FloatingTypeGPU compliance, FloatingTypeGPU dt, FloatingTypeGPU damping,
		FloatingTypeGPU* constraintLambda, FloatingTypeGPU* gradients);

	__host__  __device__ void applyToElem(PBDTetMeshFEMGPU* pTetMeshGPU, int32_t tetId, const int32_t* tetVIds,
		FloatingTypeGPU C, FloatingTypeGPU compliance, FloatingTypeGPU dt, FloatingTypeGPU damping,
		FloatingTypeGPU* constraintLambda, FloatingTypeGPU* gradients);

	//__inline__ __host__  __device__ void applyToInfiniteStiffness(TetMeshFEMGPU* pTetMeshGPU, int32_t tetId, int32_t* tetVIds,
	//	//FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
	//	FloatingTypeGPU C, FloatingTypeGPU* gradients);

	void solveMaterial(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, 
		int32_t numThreads, TetMeshFEMGPU_NeoHookean* pTetMeshGPU, cudaStream_t stream);

	void solveMaterial_oneThreadForLoop(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		int32_t numThreads, TetMeshFEMGPU_NeoHookean* pTetMeshGPU, cudaStream_t stream);

	__global__ void solveMaterialGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, int32_t* tetVIdsAll, FloatingTypeGPU* vertexInvMass,
		FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU* DmInvs, FloatingTypeGPU* devLambdas, FloatingTypeGPU dt,  
		FloatingTypeGPU devCompliance,	FloatingTypeGPU devDamping);

	__global__ void solveMaterialGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		TetMeshFEMGPU_NeoHookean* pTetMeshGPU);

	__global__ void solveNeoHookeanMaterialAggregatedGPU_kernel(TetMeshFEMGPU_NeoHookean* pTetMeshGPU);
	
	void solveNeoHookeanMaterialAggregatedGPU(int32_t numThreads, cudaStream_t stream,
		int32_t maxParallelGroupSize, TetMeshFEMGPU_NeoHookean* pTetMeshGPU);
	
	__global__ void solveMaterialGPU_oneThreadForLoop(TetMeshFEMGPU_NeoHookean tetMeshGPU);

	void solveMaterialSequantialOnCPU(TetMeshFEMGPU_NeoHookean* pTetMeshGPU);
	void solveMaterialParallelOnCPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, 
		TetMeshFEMGPU_NeoHookean* pTetMeshGPU);

#ifdef COMPILE_TEST_FUNCTIONS
	// this should only be called with one thread, one block
	void test_MaterialSolve_oneTet_call(int32_t iTet, TetMeshFEMGPU_NeoHookean* pTetMeshGPU, FloatingTypeGPU* gradDev, FloatingTypeGPU* CDev);
	__global__ void test_solveDev_oneTet(int32_t tetId, FloatingTypeGPU* vertPos, FloatingTypeGPU* vertPrevPos, const int32_t* tetVIdsAll, FloatingTypeGPU* vertexInvMass,
		FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU* DmInvs, FloatingTypeGPU* devLambdas, FloatingTypeGPU dt,
		FloatingTypeGPU devCompliance, FloatingTypeGPU devDamping, FloatingTypeGPU* gradDev, FloatingTypeGPU* CDev);

#endif // COMPILE_TEST_FUNCTIONS

	

}
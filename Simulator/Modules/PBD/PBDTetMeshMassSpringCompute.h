#pragma once
#include <CuMatrix/Buffers/ManagedBuffer.h>
#include "PBDTetMeshFEMGPU.h"

namespace GAIA {
	struct PBDTetMeshFEMGPUMassSpring : public PBDTetMeshFEMGPU
	{
		typedef std::shared_ptr<PBDTetMeshFEMGPU> SharedPtr;
		typedef PBDTetMeshFEMGPU* Ptr;

		// edge coloring
		int32_t numEdgesColoredCatergories;
		int32_t* edgesColoringEachCategorySize;
		int32_t** edgesColoringCategoriesPointers;

		FloatingTypeGPU* springLambdas;
		FloatingTypeGPU* orgLengths;
		int32_t* edges;

		// non-variant data
		FloatingTypeGPU springCompliance;
		FloatingTypeGPU springDamping;

		int32_t nEdges;

	};

	__host__ __device__ FloatingTypeGPU solveMaterialForOneEdge_MassSpring(int32_t edgeId, PBDTetMeshFEMGPUMassSpring& tetMeshGPU);

	__global__ void solveSpringConstraints_kernel(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize,
		PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);

	__global__ void solveVolumeConstraints_kernel(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);

	void solveSpringConstraintsGPU(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize,
		int32_t numThreads, PBDTetMeshFEMGPUMassSpring* pTtetMeshGPU, cudaStream_t stream);
	void solveVolumeConstraintsGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		int32_t numThreads, PBDTetMeshFEMGPUMassSpring* pTtetMeshGPU, cudaStream_t stream);

	// v2: aggregated solver version
	__global__ void solveSpringConstraintsAggregatedGPU_kernel(PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);
	void solveSpringConstraintsAggregatedGPU(int32_t numThreads, cudaStream_t stream, int32_t maxParallelGroupSize, 
		PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);

	__global__ void solveVolumeConstraintsAggregatedGPU_kernel(PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);
	void solveVolumeConstraintsAggregatedGPU(int32_t numThreads, cudaStream_t stream, int32_t maxParallelGroupSize,
		PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);

	void solveSpringConstraintsCPU(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize,
		PBDTetMeshFEMGPUMassSpring tetMeshGPU);
	void solveVolumeConstraintsCPU(const int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		PBDTetMeshFEMGPUMassSpring tetMeshGPU);

	//void solveMaterialParallelOnCPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	//	PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);

}
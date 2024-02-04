#include "PBDTetMeshMassSpringCompute.h"
#include "PBDTetMeshGeneralCompute.h"

#include "../Parallelization//GPUParallelization.h"
#include "../Parallelization//CPUParallelization.h"
#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrix.h"

#include "device_functions.h"

__inline__ __host__ __device__ FloatingTypeGPU GAIA::solveMaterialForOneEdge_MassSpring(int32_t edgeId, PBDTetMeshFEMGPUMassSpring& tetMeshGPU)
{
	//printf("Solving for edge: %d, springCompliance: %f, dt: %f\n", edgeId, tetMeshGPU.springCompliance, tetMeshGPU.dt);
	FloatingTypeGPU alpha = tetMeshGPU.springCompliance / (tetMeshGPU.dt * tetMeshGPU.dt);
	int32_t* edge = tetMeshGPU.edges + 2 * edgeId;
	int32_t v1 = edge[0];
	int32_t v2 = edge[1];
	FloatingTypeGPU* x1 = tetMeshGPU.vertPos + VERTEX_BUFFER_STRIDE * v1;
	FloatingTypeGPU* x2 = tetMeshGPU.vertPos + VERTEX_BUFFER_STRIDE * v2;

	// compute the constraint value
	FloatingTypeGPU diff[3];
	CuMatrix::vec3Minus(x1, x2, diff);
	FloatingTypeGPU l = CuMatrix::vec3Norm(diff);
	FloatingTypeGPU C = l - tetMeshGPU.orgLengths[edgeId];

	FloatingTypeGPU x1Grad[3]; 
	CuMatrix::vec3Mul(diff, 1/l, x1Grad);
	FloatingTypeGPU x2Grad[3];
	CuMatrix::vec3Mul(x1Grad, -1.f, x2Grad);

	if (l == 0.0)
	{
		return 0.f;
	}

	// apply constraints with compliance
	FloatingTypeGPU w = 0.0;
	FloatingTypeGPU damp = 0.0;

	w = tetMeshGPU.vertexInvMass[v1] + tetMeshGPU.vertexInvMass[v2];

	bool applyDamp = tetMeshGPU.springDamping != 0.f;

	if (applyDamp) {
		FloatingTypeGPU vel1[3];
		CuMatrix::vec3Minus(x1, tetMeshGPU.vertPrevPos + VERTEX_BUFFER_STRIDE * v1, vel1);
		damp += CuMatrix::vec3DotProduct(x1Grad, vel1);

		FloatingTypeGPU vel2[3];
		CuMatrix::vec3Minus(x2, tetMeshGPU.vertPrevPos + VERTEX_BUFFER_STRIDE * v2, vel2);

		damp += CuMatrix::vec3DotProduct(x2Grad, vel2);
	
		C += damp * (tetMeshGPU.springDamping) / tetMeshGPU.dt;
		w *= 1 + tetMeshGPU.springDamping / tetMeshGPU.dt;
	}

	FloatingTypeGPU dlambda = 0;
	if (w + alpha != 0) {
		dlambda = (-C - alpha * tetMeshGPU.springLambdas[edgeId]) / (w + alpha);
		//dlambda = (-C ) / (w + alpha);

	}

	if (tetMeshGPU.springCompliance != 0) {
		tetMeshGPU.springLambdas[edgeId] += dlambda;
	}

	//std::cout << "x1Grad: " << x1Grad.norm() << "\n";
	//std::cout << "x2Grad: " << x2Grad.norm() << "\n";
	//std::cout << "c: " << C << "\n";
	//std::cout << "dlambda: " << dlambda << "\n";

	/*if (isnan(x1[0]) || isnan(x2[0]) || isnan(dlambda) || isnan(invMass[v1]) || isnan(invMass[v2]) || isnan(x1Grad[0]) || isnan(x2Grad[0]))
	{
		std::cout << massSpringLambdaArr[constraintId] << "\n";
	}*/

	CuMatrix::vec3Mul(x1Grad, dlambda * tetMeshGPU.vertexInvMass[v1], x1Grad);
	CuMatrix::vec3Mul(x2Grad, dlambda * tetMeshGPU.vertexInvMass[v2], x2Grad);
	CuMatrix::vec3Add(x1, x1Grad, x1);
	CuMatrix::vec3Add(x2, x2Grad, x2);
	//x2 = x2 + (dlambda) * tetMeshGPU.vertexInvMass[v2] * (x2Grad);

	return C;
}

__global__ void GAIA::solveSpringConstraints_kernel(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	int iEdge = blockIdx.x * blockDim.x + threadIdx.x;
	if (iEdge < edgeParallelizationGroupSize)
	{
		int32_t edgeId = pEdgeParallelizationGroup[iEdge];
		solveMaterialForOneEdge_MassSpring(edgeId, *pTetMeshGPU);
	}
	return;
}

__global__ void GAIA::solveVolumeConstraints_kernel(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	int iTet = blockIdx.x * blockDim.x + threadIdx.x;
	if (iTet < tetParallelizationGroupSize)
	{
		int32_t tetId = pTetParallelizationGroup[iTet];

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

		FloatingTypeGPU* v0 = pTetMeshGPU->vertPos + id0 * 3;
		FloatingTypeGPU* v1 = pTetMeshGPU->vertPos + id1 * 3;
		FloatingTypeGPU* v2 = pTetMeshGPU->vertPos + id2 * 3;
		FloatingTypeGPU* v3 = pTetMeshGPU->vertPos + id3 * 3;

		FloatingTypeGPU gradVol[12];

		FloatingTypeGPU CVol = GAIA::solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, pTetMeshGPU->DmInvs);
		GAIA::applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, tetId, tetVIds, CVol, gradVol);
	}
	return;
}

void GAIA::solveSpringConstraintsGPU(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize,
	int32_t numThreads, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU, cudaStream_t stream)
{
	if (edgeParallelizationGroupSize > numThreads)
	{
		solveSpringConstraints_kernel KERNEL_ARGS4((edgeParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)
			(pEdgeParallelizationGroup, edgeParallelizationGroupSize, pTetMeshGPU);
	}
	else
	{
		//solveSpringConstraints_kernel KERNEL_ARGS4(1, edgeParallelizationGroupSize, 0, stream)
		//	(pEdgeParallelizationGroup, edgeParallelizationGroupSize, pTetMeshGPU);

		solveSpringConstraints_kernel KERNEL_ARGS4(1, numThreads, 0, stream)
			(pEdgeParallelizationGroup, edgeParallelizationGroupSize, pTetMeshGPU);
	}
}

void GAIA::solveVolumeConstraintsGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	int32_t numThreads, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU, cudaStream_t stream)
{
	if (tetParallelizationGroupSize > numThreads)
	{
		solveVolumeConstraints_kernel KERNEL_ARGS4((tetParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU);
	}
	else
	{
		//solveVolumeConstraints_kernel KERNEL_ARGS4(1, tetParallelizationGroupSize, 0, stream)
		//	(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU);

		solveVolumeConstraints_kernel KERNEL_ARGS4(1, numThreads, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU);
	}
}

__global__ void GAIA::solveSpringConstraintsAggregatedGPU_kernel(PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	int32_t iThread = blockIdx.x * blockDim.x + threadIdx.x;
	for (int32_t iParalGroup = 0; iParalGroup < pTetMeshGPU->numEdgesColoredCatergories; iParalGroup++)
	{
		if (iThread >= pTetMeshGPU->edgesColoringEachCategorySize[iParalGroup]) {
			return;
		}

		int32_t edgeId = pTetMeshGPU->edgesColoringCategoriesPointers[iParalGroup][iThread];

		// printf("TheadId: %d, Solving edge: %d, parallel group id %d, size: %d\n", iThread, edgeId, iParalGroup,
		// 	pTetMeshGPU->edgesColoringEachCategorySize[iParalGroup]);
		solveMaterialForOneEdge_MassSpring(edgeId, *pTetMeshGPU);
		__syncthreads();
	}
}

void GAIA::solveSpringConstraintsAggregatedGPU(int32_t numThreads, cudaStream_t stream, int32_t maxParallelGroupSize, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	solveSpringConstraintsAggregatedGPU_kernel KERNEL_ARGS4((maxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream) (pTetMeshGPU);
}

__global__ void GAIA::solveVolumeConstraintsAggregatedGPU_kernel(PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	int iThread = blockIdx.x * blockDim.x + threadIdx.x;
	for (int32_t iParalGroup = 0; iParalGroup < pTetMeshGPU->numTetsColoredCatergories; iParalGroup++)
	{
		if (iThread >= pTetMeshGPU->tetsColoringEachCategorySize[iParalGroup]) {
			return;
		}

		int32_t tetId = pTetMeshGPU->tetsColoringCategoriesPointers[iParalGroup][iThread];
		// printf("TheadId: %d, Solving tet: %d, parallel group id %d, size: %d\n", iThread, tetId, iParalGroup,
		// 	pTetMeshGPU->edgesColoringEachCategorySize[iParalGroup]);
		const int32_t* tetVIds = pTetMeshGPU->tetVIds + tetId * 4;
		int32_t id0 = tetVIds[0];
		int32_t id1 = tetVIds[1];
		int32_t id2 = tetVIds[2];
		int32_t id3 = tetVIds[3];

		FloatingTypeGPU* v0 = pTetMeshGPU->vertPos + id0 * 3;
		FloatingTypeGPU* v1 = pTetMeshGPU->vertPos + id1 * 3;
		FloatingTypeGPU* v2 = pTetMeshGPU->vertPos + id2 * 3;
		FloatingTypeGPU* v3 = pTetMeshGPU->vertPos + id3 * 3;

		FloatingTypeGPU gradVol[12];

		FloatingTypeGPU CVol = GAIA::solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, pTetMeshGPU->DmInvs);
		GAIA::applyToInfiniteStiffness(pTetMeshGPU->vertPos, pTetMeshGPU->vertexInvMass, tetId, tetVIds, CVol, gradVol);

		__syncthreads();
	}
}

void GAIA::solveVolumeConstraintsAggregatedGPU(int32_t numThreads, cudaStream_t stream, int32_t maxParallelGroupSize, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	solveVolumeConstraintsAggregatedGPU_kernel KERNEL_ARGS4((maxParallelGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream) (pTetMeshGPU);
}

void GAIA::solveSpringConstraintsCPU(int32_t* pEdgeParallelizationGroup, int32_t edgeParallelizationGroupSize, PBDTetMeshFEMGPUMassSpring tetMeshGPU)
{
	cpu_parallel_for(0, edgeParallelizationGroupSize, [&](int iEdge)
		{
			int32_t edgeId = pEdgeParallelizationGroup[iEdge];
			solveMaterialForOneEdge_MassSpring(edgeId, tetMeshGPU);
		});
}

void GAIA::solveVolumeConstraintsCPU(const int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	PBDTetMeshFEMGPUMassSpring tetMeshGPU)
{
	cpu_parallel_for(0, tetParallelizationGroupSize, [&](int iTet)
		{
			int32_t tetId = pTetParallelizationGroup[iTet];

			const int32_t* tetVIds = tetMeshGPU.tetVIds + tetId * 4;
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

			FloatingTypeGPU* v0 = tetMeshGPU.vertPos + id0 * 3;
			FloatingTypeGPU* v1 = tetMeshGPU.vertPos + id1 * 3;
			FloatingTypeGPU* v2 = tetMeshGPU.vertPos + id2 * 3;
			FloatingTypeGPU* v3 = tetMeshGPU.vertPos + id3 * 3;

			FloatingTypeGPU gradVol[12];

			FloatingTypeGPU CVol = GAIA::solveVolConstraint(tetId, v0, v1, v2, v3, gradVol, tetMeshGPU.DmInvs);
			GAIA::applyToInfiniteStiffness(tetMeshGPU.vertPos, tetMeshGPU.vertexInvMass, tetId, tetVIds, CVol, gradVol);
		});
}


#include "PBDTetMeshGeneralCompute.h"
#include "../Parallelization//GPUParallelization.h"
#include "../Parallelization//CPUParallelization.h"
#include <CuMatrix/MatrixOps/CuMatrix.h>
#include <CuMatrix/Geometry/Geometry.h>

using namespace GAIA;

using namespace CuMatrix;

void GAIA::solveBoxBoundaryConstraint(cudaStream_t stream, int32_t numVerts, int32_t numThreads, PBDTetMeshFEMGPU* pTetMeshGPU,
	FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax,
	FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic,
	FloatingTypeGPU boundaryFritionDynamic)
{

	solveBoxBoundaryConstraintGPU KERNEL_ARGS4((numVerts + numThreads - 1) / numThreads, numThreads, 0, stream) (pTetMeshGPU, xMin, xMax, yMin, yMax,
		zMin, zMax, boundaryFritionStatic, boundaryFritionDynamic);
}

void GAIA::solveBoxBoundaryConstraintOnCPU(PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax, FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic, FloatingTypeGPU boundaryFritionDynamic)
{
	FloatingTypeGPU lowerBounds[3] = { xMin, yMin, zMin };
	FloatingTypeGPU upperBounds[3] = { xMax, yMax, zMax };

	//for (size_t iVert = 0; iVert < pTetMeshGPU->nVerts; iVert++){
	cpu_parallel_for(0, pTetMeshGPU->nVerts, [&](int iVert)
		{
			solveBoxBoundaryContraintForVertex(pTetMeshGPU, pTetMeshGPU->vertPos + VERTEX_BUFFER_STRIDE * iVert, iVert, lowerBounds, upperBounds, boundaryFritionStatic, boundaryFritionDynamic);
		});
}

__host__ __device__ void GAIA::solveBoxBoundaryContraintForVertex(PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU* v,
	int32_t vId, FloatingTypeGPU* lowerBounds, FloatingTypeGPU* upperBounds, FloatingTypeGPU boundaryFritionStatic,
	FloatingTypeGPU boundaryFritionDynamic)
{
	FloatingTypeGPU vel[3];

	for (size_t iDim = 0; iDim < 3; iDim++)
	{
		FloatingTypeGPU contactNormal[3] = { 0.f, 0.f, 0.f };

		if (v[iDim] < lowerBounds[iDim])
		{
			FloatingTypeGPU penetrationDepth = lowerBounds[iDim] - v[iDim];
			v[iDim] = lowerBounds[iDim];
			
			vec3Minus(v, pTetMeshGPU->vertPrevPos + VERTEX_BUFFER_STRIDE *vId, vel);

			contactNormal[iDim] = 1;
			applyBoundaryFriction(v, penetrationDepth, contactNormal, vel, boundaryFritionStatic, boundaryFritionDynamic);

		}
		else if (v[iDim] > upperBounds[iDim])
		{
			FloatingTypeGPU penetrationDepth = v[iDim] - upperBounds[iDim];
			v[iDim] = upperBounds[iDim];

			vec3Minus(v, pTetMeshGPU->vertPrevPos + VERTEX_BUFFER_STRIDE * vId, vel);

			contactNormal[iDim] = -1;

			applyBoundaryFriction(v, penetrationDepth, contactNormal, vel, boundaryFritionStatic, boundaryFritionDynamic);

		}
	}	
}

__host__ __device__ void GAIA::applyBoundaryFriction(FloatingTypeGPU* v, FloatingTypeGPU penetrationDepth,
	FloatingTypeGPU* contactNormal, FloatingTypeGPU* vel, FloatingTypeGPU boundaryFritionStatic,
	FloatingTypeGPU boundaryFritionDynamic)
{
	// shift that is perpendicular to the contact normal
	FloatingTypeGPU diff[3]; //vel - contactNormal * (vel * contactNormal);
	FloatingTypeGPU velNormPerp = vec3DotProduct(vel, contactNormal);

	vec3Mul(contactNormal, velNormPerp, contactNormal);
	vec3Minus(vel, contactNormal, diff);

	FloatingTypeGPU dNorm = vec3Norm(diff);

	if (dNorm == 0.f)
	{
		return;
	}

	//FloatingType dynamicNorm = C * curPhysics.friction_dynamic;
	//FloatingType staticNorm = C * curPhysics.friction_static;

	FloatingTypeGPU dynamicNorm;
	if (penetrationDepth * boundaryFritionDynamic / dNorm <= 1)
	{
		dynamicNorm = penetrationDepth * boundaryFritionDynamic;
	}
	else
	{
		dynamicNorm = dNorm;
	}


	FloatingTypeGPU staticNorm = penetrationDepth * boundaryFritionStatic;

	if (dNorm > staticNorm) {
		vec3Mul(diff, dynamicNorm / dNorm, diff);
		// diff = (diff / dNorm) * dynamicNorm;
	}

	vec3Minus(v, diff, v);
}

__global__ void GAIA::solveBoxBoundaryConstraintGPU(PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax,
	FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic, FloatingTypeGPU boundaryFritionDynamic)
{
	FloatingTypeGPU lowerBounds[3] = { xMin, yMin, zMin };
	FloatingTypeGPU upperBounds[3] = { xMax, yMax, zMax };

	for (int iVert = blockIdx.x * blockDim.x + threadIdx.x;
		iVert < pTetMeshGPU->nVerts;
		iVert += blockDim.x * gridDim.x)
	{
		solveBoxBoundaryContraintForVertex(pTetMeshGPU, pTetMeshGPU->vertPos + VERTEX_BUFFER_STRIDE * iVert, iVert, lowerBounds, upperBounds,
			boundaryFritionStatic, boundaryFritionDynamic);

		// here we also prepare for the next stage: intialize the vertex inverted sign
		pTetMeshGPU->verticesInvertedSign[iVert] = false;
	}
}

__host__ __device__ void GAIA::solveInversionForOneTet(int32_t tetId, FloatingTypeGPU* vertPos, const int32_t* tetVIds, FloatingTypeGPU* vertexInvMass,
	FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU inversionSolveConstraintMultiplier)
{
	int id0 = tetVIds[0];
	int id1 = tetVIds[1];
	int id2 = tetVIds[2];
	int id3 = tetVIds[3];

	FloatingTypeGPU* v0 = vecPtr(vertPos, id0, VERTEX_BUFFER_STRIDE);
	FloatingTypeGPU* v1 = vecPtr(vertPos, id1, VERTEX_BUFFER_STRIDE);
	FloatingTypeGPU* v2 = vecPtr(vertPos, id2, VERTEX_BUFFER_STRIDE);
	FloatingTypeGPU* v3 = vecPtr(vertPos, id3, VERTEX_BUFFER_STRIDE);

	FloatingTypeGPU e1[3]; //vertexid1) - pTM->vertex(id0);
	CuMatrix::vec3Minus(v1, v0, e1);
	FloatingTypeGPU e2[3]; // = vertexid2) - pTM->vertex(id0);
	CuMatrix::vec3Minus(v2, v0, e2);
	FloatingTypeGPU e3[3]; // = vertexid3) - pTM->vertex(id0);
	CuMatrix::vec3Minus(v3, v0, e2);
	
	FloatingTypeGPU vol = CuMatrix::vec3TripleProduct(e1, e2, e3) / 6.0f; //(e1.dot(e2.cross(e3))) / 6.0f;

	if (vol <= 0.f) {
		FloatingTypeGPU grads[12];

		FloatingTypeGPU * gradC0 = grads; // = (e3 - e1).cross(e2 - e1);
		FloatingTypeGPU e3Me1[3];
		FloatingTypeGPU e2Me1[3];
		CuMatrix::vec3Minus(e3, e1, e3Me1);
		CuMatrix::vec3Minus(e2, e1, e2Me1);
		CuMatrix::vec3CrossProduct(e3Me1, e2Me1, gradC0);

		FloatingTypeGPU* gradC1 = grads + 3; // (e2.cross(e3));
		CuMatrix::vec3CrossProduct(e2, e3, gradC1);

		FloatingTypeGPU* gradC2 = grads + 6; // (e3.cross(e1));
		CuMatrix::vec3CrossProduct(e3, e1, gradC2);

		FloatingTypeGPU* gradC3 = grads + 9; // (e1.cross(e2));
		CuMatrix::vec3CrossProduct(e1, e2, gradC3);

		FloatingTypeGPU C = 6 * (vol - (1.f / tetInvRestVolume[tetId]) * inversionSolveConstraintMultiplier);
		applyToInfiniteStiffness(vertPos, vertexInvMass, tetId, tetVIds, C, grads);
	}

}

__global__ void GAIA::solveInversionConstraintGPU_kernel(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU inversionSolveConstraintMultiplier) {
	for (int iTet = blockIdx.x * blockDim.x + threadIdx.x;
		iTet < tetParallelizationGroupSize;
		iTet += blockDim.x * gridDim.x)
	{
		int32_t tetId = pTetParallelizationGroup[iTet];

		solveInversionForOneTet(tetId, pTetMeshGPU->vertPos, getTetVIdsPtr(pTetMeshGPU, tetId),
			pTetMeshGPU->vertexInvMass, pTetMeshGPU->tetInvRestVolume, inversionSolveConstraintMultiplier);
	}

}

void GAIA::solveInversionConstraintGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
	int32_t numThreads, cudaStream_t stream,
	PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU inversionSolveConstraintMultiplier)
{
	if (tetParallelizationGroupSize > numThreads)
	{
		solveInversionConstraintGPU_kernel KERNEL_ARGS4((tetParallelizationGroupSize + numThreads - 1) / numThreads, numThreads, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU, inversionSolveConstraintMultiplier);
	}
	else
	{
		solveInversionConstraintGPU_kernel KERNEL_ARGS4(1, tetParallelizationGroupSize, 0, stream)
			(pTetParallelizationGroup, tetParallelizationGroupSize, pTetMeshGPU, inversionSolveConstraintMultiplier);
	}
}

void GAIA::evaluateInversionGPU(int32_t numThreads, cudaStream_t stream, int32_t numTets, PBDTetMeshFEMGPU* pTetMeshGPU)
{
	parallel_for_tets KERNEL_ARGS4((numTets + numThreads - 1) / numThreads, numThreads, 0, stream)
		(pTetMeshGPU, [] __device__(int32_t iTet, PBDTetMeshFEMGPU * pTetMeshGPU) {
		if (CuMatrix::tetOrientedVolume(pTetMeshGPU->vertPos, getTetVIdsPtr(pTetMeshGPU, iTet)) > 0.f) {
			pTetMeshGPU->tetInvertedSign[iTet] = false;
		}
		else
		{
			pTetMeshGPU->tetInvertedSign[iTet] = true;
			// set 4 verts to inverted as well; verticesInvertedSign has been initialized previously in boundary solve
			const int32_t* tetVIds = getTetVIdsPtr(pTetMeshGPU, iTet);
			// printf("Inverted detected for tet: %d \n", iTet);
			// we are only changing from false to true, concurrency shouldn't be a problem
			pTetMeshGPU->verticesInvertedSign[tetVIds[0]] = true;
			pTetMeshGPU->verticesInvertedSign[tetVIds[1]] = true;
			pTetMeshGPU->verticesInvertedSign[tetVIds[2]] = true;
			pTetMeshGPU->verticesInvertedSign[tetVIds[3]] = true;
		}
	});
}

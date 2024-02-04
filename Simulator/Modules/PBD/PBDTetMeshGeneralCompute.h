#pragma once
#include "cuda_runtime.h"
#include <stdint.h>
#include "PBDTetMeshFEMGPU.h"
#include <CuMatrix/MatrixOps/CuMatrix.h>

namespace GAIA {
	__host__ __device__ const int32_t* getTetVIdsPtr(const int32_t* tetVIdsAll, int32_t tetId);
	__host__ __device__ const int32_t* getTetVIdsPtr(PBDTetMeshFEMGPU* pTetMeshGPU, int32_t tetId);

	void solveBoxBoundaryConstraint(cudaStream_t stream, int32_t numVerts, int32_t numThreads, PBDTetMeshFEMGPU* pTetMeshGPU,
		FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax, 
		FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic, 
		FloatingTypeGPU boundaryFritionDynamic);

	void solveBoxBoundaryConstraintOnCPU(PBDTetMeshFEMGPU* pTetMeshGPU,
		FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax,
		FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic,
		FloatingTypeGPU boundaryFritionDynamic);

	__host__ __device__ void solveBoxBoundaryContraintForVertex(PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU* v, 
		int32_t vId, FloatingTypeGPU* lowerBounds, FloatingTypeGPU* upperBounds, FloatingTypeGPU boundaryFritionStatic,
		FloatingTypeGPU boundaryFritionDynamic);
	__host__ __device__ void applyBoundaryFriction(FloatingTypeGPU* v, FloatingTypeGPU penetrationDepth,
		FloatingTypeGPU* contactNormal, FloatingTypeGPU* vel, FloatingTypeGPU boundaryFritionStatic,
		FloatingTypeGPU boundaryFritionDynamic);

	__global__ void solveBoxBoundaryConstraintGPU(PBDTetMeshFEMGPU* tetMeshGPU,
		FloatingTypeGPU xMin, FloatingTypeGPU xMax, FloatingTypeGPU yMin, FloatingTypeGPU yMax,
		FloatingTypeGPU zMin, FloatingTypeGPU zMax, FloatingTypeGPU boundaryFritionStatic,
		FloatingTypeGPU boundaryFritionDynamic);

	__host__ __device__ void solveInversionForOneTet(int32_t tetId, FloatingTypeGPU* vertPos, const int32_t* tetVIds, FloatingTypeGPU* vertexInvMass,
		FloatingTypeGPU* tetInvRestVolume, FloatingTypeGPU inversionSolveConstraintMultiplier);
	
	__host__  __device__ FloatingTypeGPU solveVolConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
		FloatingTypeGPU* gradVol, FloatingTypeGPU* DmInvs);

	__host__  __device__ void applyToInfiniteStiffness(FloatingTypeGPU* vertPos, FloatingTypeGPU* vertexInvMass, int32_t tetId, const int32_t* tetVIds,
		FloatingTypeGPU C, FloatingTypeGPU* gradients);

	__global__ void solveInversionConstraintGPU_kernel(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU inversionSolveConstraintMultiplier);

	void solveInversionConstraintGPU(int32_t* pTetParallelizationGroup, int32_t tetParallelizationGroupSize,
		int32_t numThreads, cudaStream_t stream, PBDTetMeshFEMGPU* pTetMeshGPU, FloatingTypeGPU inversionSolveConstraintMultiplier);

	void evaluateInversionGPU(int32_t numThreads, cudaStream_t stream, int32_t numTets,
		PBDTetMeshFEMGPU* pTetMeshGPU);
}

__host__ __device__ inline  FloatingTypeGPU GAIA::solveVolConstraint(int32_t tetId, FloatingTypeGPU* v0, FloatingTypeGPU* v1, FloatingTypeGPU* v2, FloatingTypeGPU* v3,
	FloatingTypeGPU* gradVol, FloatingTypeGPU* invDs)
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
	FloatingTypeGPU* invDS = invDs + tetId * 9;

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

__host__ __device__ inline void GAIA::applyToInfiniteStiffness(FloatingTypeGPU* vertPos, FloatingTypeGPU* vertexInvMass, int32_t tetId, const int32_t* tetVIds,
	FloatingTypeGPU C, FloatingTypeGPU* gradients)
{
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
	//FloatingTypeGPU damp = 0.0f;

	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		w += CuMatrix::vec3NormSquare(gs[i]) * vertexInvMass[vId];

		// FloatingTypeGPU vel[3];
		// CuMatrix::vec3Minus(pTetMeshGPU->vertPos + 3 * vId, pTetMeshGPU->vertPrevPos + 3 * vId, vel);

		// damp += CuMatrix::vec3DotProduct(gs[i], vel);
		// CPoint vel = *(mesh.verts[id]) - mesh.prevPos[id];
		// damp += gradients[i] * vel;
	}


	if (w == 0.0f) {
		return;
	}

	FloatingTypeGPU dlambda = -C / w;


	for (int32_t i = 0; i < 4; i++) {
		int32_t vId = tetVIds[i];
		FloatingTypeGPU* v = vertPos + vId * 3;

		CuMatrix::vec3Mul(gs[i], dlambda * vertexInvMass[vId], gs[i]);

		CuMatrix::vec3Add(v, gs[i], v);

		//*(mesh.verts[id]) = (*(mesh.verts[id])) + ((dlambda) * (invMass[id]) * (gradients[i]));

	}
}

__host__ __device__ inline const int32_t* GAIA::getTetVIdsPtr(const int32_t* tetVIdsAll, int32_t tetId)
{
	return tetVIdsAll + 4* tetId;
}

__host__ __device__ inline const int32_t* GAIA::getTetVIdsPtr(PBDTetMeshFEMGPU* pTetMeshGPU, int32_t tetId)
{
	return pTetMeshGPU->tetVIds + 4 * tetId;
}

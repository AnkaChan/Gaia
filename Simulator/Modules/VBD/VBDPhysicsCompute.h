#pragma once
#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrix.h"
#include "VBD_NeoHookeanGPU.h"

#define NUM_THREADS_TET_SWEEP 8
#define NUM_THREADS_VERTEX_SWEEP 32

namespace GAIA {

	void VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		cudaStream_t cudaStream);

	void VBDSolveParallelGroup_updateVertexPositionGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		int32_t numThreads, cudaStream_t cudaStream);

	void VBDSolveParallelGroup_applyAccelerationGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		int32_t numThreads, cudaStream_t cudaStream);

	void GDSolveParallelGroup_allInOneSweepGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		int32_t numThreads, cudaStream_t cudaStream);

	void GDSolveParallelGroup_BlockJacobi_allInOneSweepGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		int32_t numThreads, cudaStream_t cudaStream);

	void GDSolveParallelGroup_updatePositionSweepGPU(VBDPhysicsDataGPU* pPhysicsData, CFloatingTypeGPU stepSize, CFloatingTypeGPU acceleratorOmega, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
		int32_t numThreads, cudaStream_t cudaStream);

	void evaluateElasticEnergyGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize, FloatingTypeGPU* tetElasticEnergy,
		int32_t numThreads, cudaStream_t cudaStream);

	void VBDUpdateVelocityGPU(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices,
		int32_t numThreads, cudaStream_t cudaStream);

	__global__ void GDSolveParallelGroup_BlockJacobi_allInOneSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
		int32_t parallelGroupSize);

	__global__ void GDSolveParallelGroup_updatePositionSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, CFloatingTypeGPU stepSize, CFloatingTypeGPU acceleratorOmega, int32_t* parallelGroupsHead,
		int32_t parallelGroupSize);

	// this is the recommended version, which does the block-thread two layer parallelism
	__global__ void VBDSolveParallelGroup_allInOne_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
		int32_t parallelGroupSize);
	// copy the updated position to the original position buffer to ensure GS iteration
	// also in charge of acceleration
	__global__ void VBDSolveParallelGroup_updateVertexPosition_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
		int32_t parallelGroupSize);
	__global__ void VBDSolveParallelGroup_applyAcceleration_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
		int32_t parallelGroupSize);

	__global__ void evaluateElasticEnergy_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize, FloatingTypeGPU* tetElasticEnergy);

	__global__ void VBDUpdateVelocity_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices);
	
	__global__ void VBDApplyAccelerator_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices, FloatingTypeGPU acceleratorOmega);

	__global__ void GDSolveParallelGroup_vertexSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	int32_t parallelGroupSize);

	// compute and accumulate the force and hessian, used by all in one sweep
	__device__ __host__  void accumulateMaterialForceAndHessianForVertex_NeoHookean(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t iV,
		FloatingTypeGPU* force, FloatingTypeGPU* h, CFloatingTypeGPU dt);

	GPU_CPU_INLINE_FUNC  void caculateMaterialForceAndHessian_NeoHookean_oneTet(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU,
		int32_t tetId, int32_t vertexOrder, FloatingTypeGPU* force, FloatingTypeGPU* h, CFloatingTypeGPU dt);

	// compute the force and hessian on the fly
	GPU_CPU_INLINE_FUNC void accumlateCollisionForceAndHessianPerCollision_v2(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
		const CollisionDataGPU& collisionResult, int32_t collisionResultId, int32_t collisionVertexOrder, FloatingTypeGPU* force, FloatingTypeGPU* hessian);

	GPU_CPU_INLINE_FUNC  void dampVelocity(FloatingTypeGPU* velocity, CFloatingTypeGPU maxVelocityMagnitude,
		CFloatingTypeGPU exponentialVelDamping, CFloatingTypeGPU constantVelDamping);

	GPU_CPU_INLINE_FUNC  void updateVertexVelocity(VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU dt);
	
	GPU_CPU_INLINE_FUNC  void evaluateNeoHookeanMaterialForceAndHessian(CFloatingTypeGPU miu, CFloatingTypeGPU lmbd, CFloatingTypeGPU A, int32_t tetId, int32_t vOrderInTet,
		const int32_t* tetVIds, CFloatingTypeGPU* vs, CFloatingTypeGPU* vsPrev, CFloatingTypeGPU* DmInvs, FloatingTypeGPU* force, FloatingTypeGPU* h,
		CFloatingTypeGPU dampingVolume, CFloatingTypeGPU dampingShear, CFloatingTypeGPU dt);
	
	GPU_CPU_INLINE_FUNC  void evaluateNeoHookeanMaterialForceAndDiagonalHessian(FloatingTypeGPU miu, FloatingTypeGPU lmbd, FloatingTypeGPU A, int32_t tetId, int32_t vOrderInTet,
		const int32_t* tetVIds, FloatingTypeGPU* vs, CFloatingTypeGPU* DmInvs, FloatingTypeGPU* force, FloatingTypeGPU* h);
	
	GPU_CPU_INLINE_FUNC  void accumulateCollisionForceAndHessian(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU, int32_t iV,
	FloatingTypeGPU* force, FloatingTypeGPU* h);

	GPU_CPU_INLINE_FUNC  CFloatingTypeGPU evaluateNeoHookeanEnergy(FloatingTypeGPU miu, FloatingTypeGPU lmbd, FloatingTypeGPU A, int32_t tetId, 
	const int32_t* tetVIds, FloatingTypeGPU* vs, CFloatingTypeGPU* DmInvs);

	GPU_CPU_INLINE_FUNC void GDStep_vertexSweep_allInOne(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId);
	
	GPU_CPU_INLINE_FUNC  void GDStep_vertexSweep_allInOne_blockJacobi(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId);


	//void VBDSolveParallelGroup_tetSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead, int32_t tetParallelGroupSize,
	//	int32_t numThreads, cudaStream_t cudaStream);
	//
	//// update and precompute the collision info
	//void VBDSolveParallelGroup_collisionSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t numActiveCollisions, int32_t* activeCollisionsEachParallelGroup,
	//	int32_t numThreads, cudaStream_t cudaStream);

	//void VBDSolveParallelGroup_vertexSweepAcceleratedGS(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize,
	//	FloatingTypeGPU acceleratorOmega, int32_t numThreads, cudaStream_t cudaStream);
	
	//void VBDSolveParallelGroup_allInOneSweep(VBDPhysicsDataGPU* pPhysicsData, int32_t* vertexParallelGroupsHead, int32_t veretxParallelGroupSize, FloatingTypeGPU acceleratorOmega,
	//	int32_t numThreads, cudaStream_t cudaStream);

	//__global__ void VBDSolveParallelGroup_allInOne_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead, 
	//	int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega);

	//__global__ void VBDSolveParallelGroup_collisionSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t numActiveCollisions, 
	//	int32_t* activeCollisionsEachParallelGroup);

	//__global__ void VBDSolveParallelGroup_tetSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead,
	//	int32_t tetParallelGroupSize);

	//__global__ void VBDSolveParallelGroup_tetSweep_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead,
	//	int32_t tetParallelGroupSize);

	//__global__ void GDSolveParallelGroup_tetSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* tetParallelGroupsHead,
	//	int32_t tetParallelGroupSize);

	//__global__ void VBDSolveParallelGroup_vertexSweep_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	//	int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega = 1.0f);

	//__global__ void VBDSolveParallelGroup_vertexSweep_kernel_V2(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	//	int32_t parallelGroupSize, CFloatingTypeGPU acceleratorOmega = 1.0f);

	//__global__ void VBDSolveParallelGroup_vertexSweepAcceleratedGS_kernel(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	//	int32_t parallelGroupSize, FloatingTypeGPU acceleratorOmega);



	//__global__ void VBDSolveParallelGroup_kernel_serial(VBDPhysicsDataGPU* pPhysicsData, int32_t* parallelGroupsHead,
	//	int32_t parallelGroupSize);

	//__global__ void VBDUpdateRelativeVelocity_kernel(VBDPhysicsDataGPU* pPhysicsData, const int32_t* vertexAllParallelGroups, const int32_t numVertices);

	//// accumulate the force and hessian computed from the previous tet sweep, used by tet-vertex two stage sweep
	//GPU_CPU_INLINE_FUNC  void accumulateMaterialForceAndHessianForVertex_NeoHookean_preCompuated(const VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t iV,
	//	FloatingTypeGPU* force, FloatingTypeGPU* h);

	//__device__ __host__ void updateCollisionInfoGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t meshId, int32_t vertexId);
	//__device__ __host__ void computeRelativeVelocityGPU(VBDPhysicsDataGPU* pPhysicsData, int32_t meshId, int32_t vertexId);

	//// pTetMeshGPU: the mesh that collisionResult belongs to
	//// collisionResultId: the primitive id of the collisionResult
	//// collisionVertexOrder: 0~3
	//GPU_CPU_INLINE_FUNC void accumlateCollisionForceAndHessianPerCollision(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
	//	const CollisionDataGPU& collisionResult, int32_t collisionResultId, int32_t collisionVertexOrder, FloatingTypeGPU* force, FloatingTypeGPU* hessian);

	//GPU_CPU_INLINE_FUNC  FloatingTypeGPU computeCollisionForcePerCollision(const VBDPhysicsDataGPU* pPhysicsData, const VBDBaseTetMeshGPU* pTetMeshGPU,
	//	const CollisionDataGPU& collisionResult, int32_t collisionResultId, int32_t collisionVertexOrder,  FloatingTypeGPU* contactNormal);

	//GPU_CPU_INLINE_FUNC  void VBDStep_allInOne(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId);
	//__inline__ __device__   void VBDStep_vertexSweep(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU acceleratorOmega = 1.0f);

	////GPU_CPU_INLINE_FUNC  void VBDStep_vertexSweep(VBDPhysicsDataGPU* pPhysics, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, int32_t vertexId, CFloatingTypeGPU acceleratorOmega = 1.0f);

	//GPU_CPU_INLINE_FUNC  void applyBoundaryFriction(VBDPhysicsDataGPU* pPhysics, VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId,
	//	CFloatingTypeGPU dt, CFloatingTypeGPU vertInvMass);

	//GPU_CPU_INLINE_FUNC  void applyCollisionFriction(VBDPhysicsDataGPU* pPhysics, VBDBaseTetMeshGPU* pTetMeshGPU, int32_t vertexId,
	//	CFloatingTypeGPU dt, CFloatingTypeGPU vertInvMass);

	//GPU_CPU_INLINE_FUNC  void applyFrictionalVelocityDampingBoundary(FloatingTypeGPU* velocity, CFloatingTypeGPU* contactNormal, CFloatingTypeGPU vertInvMass,
	//	CFloatingTypeGPU frictionRatio, CFloatingTypeGPU contactForce, CFloatingTypeGPU dt);

	//GPU_CPU_INLINE_FUNC  void applyFrictionalVelocityDamping(FloatingTypeGPU* velocity, CFloatingTypeGPU* tangentialVelDirection,
	//	CFloatingTypeGPU tangentialVelNorm, CFloatingTypeGPU contactForceNorm, CFloatingTypeGPU vertInvMass, CFloatingTypeGPU frictionRatio,  CFloatingTypeGPU dt);



}
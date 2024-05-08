#pragma once

#include "../TetMesh/TetMeshFEMGPU.h"

#define GPU_COLLISION_RELATION_PREALLOCATION_SIZE 4
#define GPU_JACOBI_DX

namespace GAIA {
	struct CollisionRelationGPU
	{
		int32_t meshId;
		int32_t collisionPrimitiveId;   // eithe vertex id or edge id
		//IdType collisionId;	  // must be 0
		int32_t collisionType;     // 0: v-f; 1: edge edge
		int32_t collisionVertexOrder; // 0~2 for v-f contact, 0~3 for e-e contact
	};

	struct CollisionDataGPU
	{
		CollisionRelationGPU collisionRelations[GPU_COLLISION_RELATION_PREALLOCATION_SIZE];
		// outputs
		//FloatingTypeGPU closestPointNormal[3];
		//FloatingTypeGPU closestSurfacePt[3];
		FloatingTypeGPU closestSurfacePtBarycentrics[3];

		// when computing collision infos: 
		//     0~2: force
		//     3~11: hessian
		// when applying frictional vel damping: 
		//     0~2: tangential relative velocity direction
		//     3: tangential relative velocity norm
		//     4: tangential relative force norm
		FloatingTypeGPU collisionForceAndHessian[12]; 
		// inputs
		int32_t closestSurfaceFaceVIds[3];
		int32_t closestSurfaceFaceId;
		int32_t intersectedTMeshId;
		// 0 or 1, records whether the point is the v side of a v-f collision
		// 0 doesn't mean it's not colliding, it mean it's not from the v side of a v-f collision
		int32_t activeColliding;       
		int32_t numCollisionRelations; // < GPU_COLLISION_RELATION_PREALLOCATION_SIZE 

	};

	struct VBDBaseTetMeshGPU : TetMeshFEMGPU {
		// 12 * numTets: 0~2 force, 3~11 hessian
		FloatingTypeGPU* tetForceAndHessians;
#ifdef GPU_JACOBI_DX
		// for jacobi style iteration
		FloatingTypeGPU* positionsNew;
#endif // !GPU_JACOBI_DX

		// for Chebyshev accelerator
		FloatingTypeGPU* positionsPrevIter;
		int8_t* activeCollisionMask;
		FloatingTypeGPU* inertia;

		// numSurfaceVertices
		CollisionDataGPU* collisionData;
		// GPU_COLLISION_RELATION_PREALLOCATION_SIZE * numVertices

		FloatingTypeGPU maxVelocityMagnitude;
		FloatingTypeGPU exponentialVelDamping;
		FloatingTypeGPU constantVelDamping;
		FloatingTypeGPU frictionDynamic;
		FloatingTypeGPU frictionEpsV;

		GPU_CPU_INLINE_FUNC FloatingTypeGPU* getVert(int32_t iVert) { return vertPos + VERTEX_BUFFER_STRIDE * iVert; }
		GPU_CPU_INLINE_FUNC CFloatingTypeGPU* getVert(int32_t iVert) const { return vertPos + VERTEX_BUFFER_STRIDE * iVert; }

		GPU_CPU_INLINE_FUNC FloatingTypeGPU* getVertPrevPos(int32_t iVert) { return vertPrevPos + VERTEX_BUFFER_STRIDE * iVert; }
		GPU_CPU_INLINE_FUNC CFloatingTypeGPU* getVertPrevPos(int32_t iVert) const { return vertPrevPos + VERTEX_BUFFER_STRIDE * iVert; }

		GPU_CPU_INLINE_FUNC CollisionDataGPU& getCollisionData(int32_t iVert) { return collisionData[iVert]; }
		GPU_CPU_INLINE_FUNC const CollisionDataGPU& getCollisionData(int32_t iVert) const { return collisionData[iVert]; }
		
		GPU_CPU_INLINE_FUNC CollisionRelationGPU& getCollisionRelation(int32_t iVert, int32_t iRelation) 
		{ 
			return collisionData[iVert].collisionRelations[iRelation];
		}

		GPU_CPU_INLINE_FUNC const CollisionRelationGPU& getCollisionRelation(int32_t iVert, int32_t iRelation) const
		{
			return collisionData[iVert].collisionRelations[iRelation];
		}

		GPU_CPU_INLINE_FUNC int32_t numCollisionRelation(int32_t iVert) const { return collisionData[iVert].numCollisionRelations; }
	};

	struct VBDPhysicsDataGPU : PhysicsDataGPU
	{
		VBDBaseTetMeshGPU** tetMeshes;
		// non-variant data
		int32_t numMeshes;
		FloatingTypeGPU worldBounds[6];
		// int32_t** vertexParallelGroups; // can be passed through the kernel call
		FloatingTypeGPU boundaryCollisionStiffness;
		FloatingTypeGPU boundaryFrictionDynamic;
		FloatingTypeGPU boundaryFrictionEpsV;

		FloatingTypeGPU collisionStiffness;
		FloatingTypeGPU collisionAirDistance;
		FloatingTypeGPU dtSqrReciprocal;
		FloatingTypeGPU stepSize;
		FloatingTypeGPU stepSizeGD;
		//FloatingTypeGPU solveOffHeight;
		
		int32_t useAccelerator;
		int32_t useBlockJacobi; // for GD only
	};

}

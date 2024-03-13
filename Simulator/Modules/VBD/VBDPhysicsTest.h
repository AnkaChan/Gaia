#pragma once

#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrix.h"
#include "VBD_BaseMeshGPU.h"

namespace GAIA {
	void printGPUCollisionData(int32_t numActiveCollisions, int32_t* activeCollisionsEachParallelGroup, VBDPhysicsDataGPU* pPhysics);

	GPU_CPU_INLINE_FUNC void printGPUCollisionData(const CollisionDataGPU* collisionData) {
		printf("    collide with (%d, %d, %d) from mesh %d:\n",
			collisionData->closestSurfaceFaceVIds[0], collisionData->closestSurfaceFaceVIds[1], collisionData->closestSurfaceFaceVIds[2],
			collisionData->intersectedTMeshId);

		//printf("    closestSurfacePt (%f, %f, %f)\n", collisionData->closestSurfacePt[0],
		//	collisionData->closestSurfacePt[1], collisionData->closestSurfacePt[2]);
		//printf("    closestPointNormal (%f, %f, %f)\n", collisionData->closestPointNormal[0],
		//	collisionData->closestPointNormal[1], collisionData->closestPointNormal[2]);
		printf("    closestSurfacePtBarycentrics (%f, %f, %f)\n", collisionData->closestSurfacePtBarycentrics[0],
			collisionData->closestSurfacePtBarycentrics[1], collisionData->closestSurfacePtBarycentrics[2]);
	}

	GPU_CPU_INLINE_FUNC void printVBDPhysicsDataGPUForVertex(const VBDPhysicsDataGPU* pPhysicsData, 
		int32_t iMesh, int32_t iV) {
		const VBDBaseTetMeshGPU* pTetMeshGPU = pPhysicsData->tetMeshes[iMesh];
		const CollisionDataGPU& collisionData = pTetMeshGPU->getCollisionData(iV);
		printf("------Print collision info for vertex %d of mesh %d------\n", iMesh, iV);

		if (collisionData.activeColliding)
		{
			//printf("collision handling for vertex: %d, collisionOrder: %d\n", iV, 3);
			//printf("    collision info:\n");
			//printf("    collide with (%d, %d, %d) from mesh %d:\n",
			//	collisionData.closestSurfaceFaceVIds[0], collisionData.closestSurfaceFaceVIds[1], collisionData.closestSurfaceFaceVIds[2],
			//	collisionData.intersectedTMeshId);
			printf("v side:\n");
			printGPUCollisionData(&collisionData);

		}

		if (collisionData.numCollisionRelations)
		{
			for (int32_t iCollision = 0; iCollision < collisionData.numCollisionRelations; iCollision++)
			{
				printf("f side (%d of %d):\n", iCollision, collisionData.numCollisionRelations);

				const CollisionRelationGPU& collisionRelation = collisionData.collisionRelations[iCollision];
				const VBDBaseTetMeshGPU* pTetMeshGPUOther = pPhysicsData->tetMeshes[collisionRelation.meshId];
				const CollisionDataGPU& collisionDataOther = pTetMeshGPUOther->getCollisionData(collisionRelation.collisionPrimitiveId);

				/*	printf("collision handling for vertex: %d, collisionPrimitiveId: d, collisionOrder: %d\n",
						iV, collisionRelation.collisionPrimitiveId, collisionRelation.collisionVertexOrder);
					printf("    collision info:\n");
					printf("    collide with (%d, %d, %d) from mesh %d:\n",
						collisionDataOther.closestSurfaceFaceVIds[0], collisionDataOther.closestSurfaceFaceVIds[1], collisionDataOther.closestSurfaceFaceVIds[2],
						collisionDataOther.intersectedTMeshId);*/

				printGPUCollisionData(&collisionDataOther);
			}
		}
	}



}
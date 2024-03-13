#pragma once
#include "../Types/Types.h"
#include "VBD_BaseMaterial.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

#define COLLISION_RELATION_PREALLOCATE 4

namespace GAIA {
	struct CollisionRelation
	{
		IdType meshId;
		IdType surfaceVertexId;          // surface VId of surface mesh
		IdType collisionId;
		IdType collisionType;     // 0: v-f; 1: edge edge
		IdType collisionVertexOrder; // 0 ~ 2 for v-f contact, 0~4 for e-e contact
	};

	typedef CPArray<CollisionRelation, COLLISION_RELATION_PREALLOCATE> CollisionRelationList;

	struct ActiveCollisionList
	{
		void initialize(const std::vector<std::shared_ptr<VBDBaseTetMesh>>& basetetMeshes,
			const std::vector<std::vector<IdType>>& vertexParallelGroups, CFloatingType activeCollisionListPreAllocationRatio)
		{
			tMeshes = basetetMeshes;
			vertexCollisionRelations.resize(basetetMeshes.size());
			for (size_t iMesh = 0; iMesh < basetetMeshes.size(); iMesh++)
			{
				vertexCollisionRelations[iMesh].resize(basetetMeshes[iMesh]->numVertices());
			}

			initializeParallelGroup(vertexParallelGroups, activeCollisionListPreAllocationRatio);
		}

		void initializeParallelGroup(const std::vector<std::vector<IdType>>& vertexParallelGroups, CFloatingType activeCollisionListPreAllocationRatio);
		

		void clear();

		void addToActiveCollisionList(CollisionDetectionResult& collisionResult, size_t numIntersections = -1);

		std::vector<std::shared_ptr<VBDBaseTetMesh>> tMeshes;
		// nCollision x (iMesh, vertexTetMeshId)
		std::vector<std::pair<IdType, IdType>> activeCollisions;
		// records how each vertex is connected to other vertices because of collision
		std::vector<std::vector<CollisionRelationList>> vertexCollisionRelations;

		// for parallel collision handling
		// 2 * nVerticesEachGroup; meshid1, vertexId1, meshId2, vertexId2,
		// reserved the size of each group to make sure it's always enough
		std::vector<VecDynamicI> activeCollisionsEachParallelGroup;
		std::vector<size_t> numActiveCollisionsEachParallelGroup;
	};

}
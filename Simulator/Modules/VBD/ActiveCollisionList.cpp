#include "ActiveCollisionList.h"
#include "../Parallelization/CPUParallelization.h"

using namespace GAIA;

void GAIA::ActiveCollisionList::initializeParallelGroup(const std::vector<std::vector<IdType>>& vertexParallelGroups, CFloatingType activeCollisionListPreAllocationRatio)
{
	size_t numberOfParallelGroups = vertexParallelGroups.size();
	activeCollisionsEachParallelGroup.resize(numberOfParallelGroups);
	numActiveCollisionsEachParallelGroup.resize(numberOfParallelGroups);
	for (int iGroup = 0; iGroup < numberOfParallelGroups; iGroup++)
	{
		auto& currentVertexColorGroup = vertexParallelGroups[iGroup];
		size_t activeCollisionsEachParallelGroupSize = size_t(currentVertexColorGroup.size() * activeCollisionListPreAllocationRatio);
		if (activeCollisionsEachParallelGroupSize % 2)
		{
			activeCollisionsEachParallelGroupSize++;
		}
		activeCollisionsEachParallelGroup[iGroup].resize(activeCollisionsEachParallelGroupSize);

	}
}

void GAIA::ActiveCollisionList::clear()
{
	activeCollisions.clear();
	cpu_parallel_for(0, vertexCollisionRelations.size(), [&](int iMesh) {
		std::vector<CollisionRelationList>& meshColRelationList = vertexCollisionRelations[iMesh];
		for (size_t iVert = 0; iVert < meshColRelationList.size(); iVert++)
		{
			vertexCollisionRelations[iMesh][iVert].clear();
		}

		});

	for (size_t iGroup = 0; iGroup < numActiveCollisionsEachParallelGroup.size(); iGroup++)
	{
		numActiveCollisionsEachParallelGroup[iGroup] = 0;
	}

}

void GAIA::ActiveCollisionList::addToActiveCollisionList(CollisionDetectionResult& collisionResult, size_t numIntersections)
{
	if (numIntersections == -1)
	{
		numIntersections = collisionResult.numIntersections();
	}
	int iMesh = collisionResult.idTMQuery;
	int surfaceVIdTetMesh = collisionResult.idVQuery;
	bool isValidCollision = false;
	for (size_t iIntersection = 0; iIntersection < numIntersections; iIntersection++)
	{
		CollidingPointInfo& collidingPt = collisionResult.collidingPts[iIntersection];
		if (collidingPt.shortestPathFound)
		{
			isValidCollision = true;
			CPArray<int, 4> parallelGroupsAdded;
			// add the collision Id to the active collision list of the vertex' parallel group

			VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
			int iParallelGroupV = pTetMesh->vertexParallelGroups[surfaceVIdTetMesh];
			int activeColId = numActiveCollisionsEachParallelGroup[iParallelGroupV]; // can be replaced by atom operation
			numActiveCollisionsEachParallelGroup[iParallelGroupV]++;				 // can be replaced by atom operation
			if (activeCollisionsEachParallelGroup[iParallelGroupV].size() > (activeColId * 2))
			{
				activeCollisionsEachParallelGroup[iParallelGroupV](activeColId * 2) = iMesh;
				activeCollisionsEachParallelGroup[iParallelGroupV](activeColId * 2 + 1) = surfaceVIdTetMesh;
				parallelGroupsAdded.push_back(iParallelGroupV);
			}
			else
			{
				std::cout << "Parallel group: " << iParallelGroupV << " has exceeded the max number collisions! Skipping more collisions\n";
			}

			const int meshId_intersecting = collidingPt.intersectedMeshId;
			const int surfaceFaceId = collidingPt.closestSurfaceFaceId;
			VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();

			for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
			{
				IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);

				// add collision Id to the active collision list of the  the face vertex's group
				int iParallelGroupFV = pTetMesh_intersecting->vertexParallelGroups[faceVId];
				// avoid adding to the same group more than once
				if (!parallelGroupsAdded.has(iParallelGroupFV))
				{
					int activeColIdFV = numActiveCollisionsEachParallelGroup[iParallelGroupFV];
					if (activeCollisionsEachParallelGroup[iParallelGroupFV].size() > (activeColIdFV * 2))
					{
						numActiveCollisionsEachParallelGroup[iParallelGroupFV]++;
						activeCollisionsEachParallelGroup[iParallelGroupFV](activeColIdFV * 2) = iMesh;
						activeCollisionsEachParallelGroup[iParallelGroupFV](activeColIdFV * 2 + 1) = surfaceVIdTetMesh;
						parallelGroupsAdded.push_back(iParallelGroupFV);
					}
					else
					{
						std::cout << "Parallel group: " << iParallelGroupFV << " has exceeded the max number collisions! Skipping more collisions\n";
					}

				}
			}
		}
	}
	if (isValidCollision)
	{
		activeCollisions.push_back({ iMesh, surfaceVIdTetMesh });
	}
}
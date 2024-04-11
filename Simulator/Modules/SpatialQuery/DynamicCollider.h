#pragma once
#include "MeshClosestPointQuery.h"
#include "ColiiderTriMeshBase.h"
namespace GAIA {
	struct DynamicColliderParameters : public MeshClosestPointQueryParameters
	{
		typedef std::shared_ptr<DynamicColliderParameters> SharedPtr;
		typedef DynamicColliderParameters* Ptr;

		DynamicColliderParameters()
		{
			maxQueryDis = 0.5f;
		}


	};


	struct DynamicCollider
	{
		typedef std::shared_ptr<DynamicCollider> SharedPtr;
		typedef DynamicCollider* Ptr;

		DynamicCollider(const DynamicColliderParameters::SharedPtr in_pQueryParams)
			: pQueryParams(in_pQueryParams)
		{
		};

		void initialize(std::vector<ColliderTrimeshBase::SharedPtr>& meshes_in)
		{
			colliderMeshes = meshes_in;
		}

		void updateColliderMeshes(IdType frameId, IdType substep, IdType iter)
		{
			for (int iColliderMesh = 0; iColliderMesh < colliderMeshes.size(); iColliderMesh++)
			{
				colliderMeshes[iColliderMesh]->update(frameId, substep, iter);
			}
		}

		const DynamicColliderParameters::SharedPtr pQueryParams;
		const BasePhysicsParams::SharedPtr pPhysicsParams;
		std::vector<ColliderTrimeshBase::SharedPtr> colliderMeshes;
	};

}
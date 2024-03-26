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

		void updateColliderMeshes()
		{
		}

		const DynamicColliderParameters::SharedPtr pQueryParams;
		std::vector<ColliderTrimeshBase::SharedPtr> colliderMeshes;
	};

}
#pragma once
#include "MeshClosestPointQuery.h"



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
			mQuery = std::make_shared<MeshClosestPointQuery>(pQueryParams);
		};

		bool collideQuery(const Vec3 p, TriMeshClosestPointQueryResult* pClosestPtResult, bool computeNormal = false)
		{
			return mQuery->closestPointQuery(p, pClosestPtResult, computeNormal);
		}

		void initialize(std::vector<TriMeshFEM::SharedPtr> meshes)
		{
			return mQuery->initialize(meshes);
		}

		void updateBVH()
		{
			return mQuery->updateBVH();
		}

		const DynamicColliderParameters::SharedPtr pQueryParams;
		MeshClosestPointQuery::SharedPtr mQuery;
	};

}
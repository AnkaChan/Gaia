#pragma once

#include <embree3/rtcore.h>

#include "../TriMesh/TriMesh.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

#define NUM_QUERY_PRE_ALLOC 2

namespace GAIA {
	struct MeshClosestPointQuery;

	struct ClosestPointInfo
	{
		int queryPointId = -1;
		int closestFaceId = -1;
		int closestMeshId = -1;
		ClosestPointOnTriangleType closestPtType;
		ClosestPointOnPrimitiveType primitiveType;
		// the id of the primitive the cloesest point lies on
		// depends on closestPtType, it can be a vertex, an edge or a face
		int primitiveId = -1;
		FloatingType d = 0.f;
		Vec3 closestPt;
		Vec3 closestPtNormal;
		Vec3 barycentrics;
	};

	struct TriMeshClosestPointQueryResult {
		// configs
		bool found = false;
		bool computeNormal = false;
		MeshClosestPointQuery* pMatcher = nullptr;

		// outputs
		CPArray<ClosestPointInfo, NUM_QUERY_PRE_ALLOC> contactPts;

		void reset() {
			found = false;
			contactPts.clear();
		}
	};


	struct MeshClosestPointQueryParameters : MF::BaseJsonConfig
	{
		typedef std::shared_ptr<MeshClosestPointQueryParameters> SharedPtr;
		typedef MeshClosestPointQueryParameters* Ptr;

		bool fromJson(nlohmann::json& j) override;
		bool toJson(nlohmann::json& j) override;

		FloatingType maxQueryDis = 1.2;
		// if set to false it will record all the points inside the search radius
		// internal parameter, cannot be set from json
		bool onlyKeepClosest = true; 
	};

	struct MeshClosestPointQuery
	{
		typedef std::shared_ptr<MeshClosestPointQuery> SharedPtr;
		typedef MeshClosestPointQuery* Ptr;

		MeshClosestPointQuery(const MeshClosestPointQueryParameters::SharedPtr in_pParams);

		~MeshClosestPointQuery()
		{
			rtcReleaseScene(targetMeshFacesScene);
		}

		void updateBVH(RTCBuildQuality tetMeshSceneQuality = RTC_BUILD_QUALITY_REFIT);
		void initialize(std::vector<TriMeshFEM::SharedPtr>  in_targetMeshes);
		bool closestPointQuery(const Vec3 p, TriMeshClosestPointQueryResult* pClosestPtResult, bool computeNormal = false);


	public:
		const MeshClosestPointQueryParameters::SharedPtr pParams;
		std::vector<TriMeshFEM::SharedPtr> targetMeshes;

		RTCScene targetMeshFacesScene;
		RTCDevice device;
	};
}
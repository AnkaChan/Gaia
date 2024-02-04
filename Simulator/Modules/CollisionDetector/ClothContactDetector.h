#pragma once
#include "../SpatialQuery/MeshClosestPointQuery.h"

namespace GAIA {
	struct ClothContactDetector;
	struct ClothContactDetectorParameters : public MeshClosestPointQueryParameters
	{
		typedef std::shared_ptr<ClothContactDetectorParameters> SharedPtr;
		typedef ClothContactDetectorParameters* Ptr;
		ClothContactDetectorParameters() {
			// we need to record all the points inside the seach radius to caculate the contact force
			onlyKeepClosest = false;
			maxQueryDis = 0.4f;
		}

		bool fromJson(nlohmann::json& j) override;
		bool toJson(nlohmann::json& j) override;
	};

	struct ClothVFContactQueryResult : public TriMeshClosestPointQueryResult {
		// query point info
		TriMeshFEM* pMeshQuery = nullptr;
		int queryPrimitiveId = -1;

	};

	struct EEClosestPointInfo
	{
		int closestEdgeId = -1;
		int closestMeshId = -1;

		FloatingType mu_this = 0.f;
		FloatingType mu_opposite = 0.f;

		Vec3 c1, c2; // c1 is from the current edge and c2 is from the opposite edge

		FloatingType d = 0.f;
	};

	struct ClothEEContactQueryResult {
		CPArray<EEClosestPointInfo, NUM_QUERY_PRE_ALLOC> contactPts;

		// query point info
		TriMeshFEM* pMeshQuery = nullptr;
		ClothContactDetector* pContactDetector = nullptr;

		int queryPrimitiveId = -1;
		bool found = true;

		void reset() {
			found = false;
			contactPts.clear();
		}
	};

	struct ClothContactDetector : public MeshClosestPointQuery
	{
		typedef std::shared_ptr<ClothContactDetector> SharedPtr;
		typedef ClothContactDetector* Ptr;

		ClothContactDetector(const ClothContactDetectorParameters::SharedPtr pParameters);
		void initialize(std::vector<TriMeshFEM::SharedPtr> in_targetMeshes);
		bool contactQueryVF(TriMeshFEM* pClothMesh, IdType vId, ClothVFContactQueryResult* pResult);
		bool contactQueryEE(TriMeshFEM* pClothMesh, IdType vId, ClothEEContactQueryResult* pResult);

		void updateBVH(RTCBuildQuality sceneQuality = RTC_BUILD_QUALITY_REFIT);


		const ClothContactDetectorParameters& parameters() const
		{
			return *(ClothContactDetectorParameters*)pParams.get();
		}

		RTCScene targetMeshEdgesScene;

	};
}
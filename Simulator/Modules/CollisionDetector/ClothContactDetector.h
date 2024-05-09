#pragma once

#define FV_CONTACT_PREALLOCATE 32
#define VF_CONTACT_PREALLOCATE 8
#define EE_CONTACT_PREALLOCATE 8
#include <embree3/rtcore.h>

#include "../TriMesh/TriMesh.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

namespace GAIA {
	struct ClothContactDetector;
	struct ClothContactDetectorParameters : public MF::BaseJsonConfig
	{
		typedef std::shared_ptr<ClothContactDetectorParameters> SharedPtr;
		typedef ClothContactDetectorParameters* Ptr;
		ClothContactDetectorParameters() {
			// we need to record all the points inside the seach radius to caculate the contact force
			maxQueryDis = 0.4f;
		}

		// if true, the contact detector will construct a BVH for the vertices
		// FV query is equivalent to VF query, but it can be more efficient in some cases
		bool supportFVQuery = false;
		FloatingType maxQueryDis = 1.2;
		// bool caculateFaceMinDis = false;

		bool fromJson(nlohmann::json& j) override;
		bool toJson(nlohmann::json& j) override;
	};

	struct VFContactPointInfo
	{
		int contactVertexId = -1;
		int contactVertexSideMeshId = -1;
		int contactFaceId = -1;
		int contactFaceSideMeshId = -1;
		ClosestPointOnTriangleType closestPtType;
		ClosestPointOnPrimitiveType primitiveType;
		// the id of the primitive the cloesest point lies on
		// depends on closestPtType, it can be a vertex, an edge or a face
		int primitiveId = -1;
		FloatingType d = 0.f;
		Vec3 contactPoint;
		Vec3 contactPointNormal;
		Vec3 barycentrics;

	};

	struct ClothVFContactQueryResult {
		// configs
		bool found = false;
		bool computeNormal = false;
		ClothContactDetector* pContactDetector = nullptr;

		ClothVFContactQueryResult() {
			minDisToPrimitives = std::numeric_limits<FloatingType>::max();
		}

		// VF contact can be detected from both VF and FV query
		int queryPrimitiveId = -1; // face id for FV query and vertex id for VF query
		int queryMeshId = -1;

		// outputs
		CPArray<VFContactPointInfo, VF_CONTACT_PREALLOCATE> contactPts;

		void reset() {
			minDisToPrimitives = std::numeric_limits<FloatingType>::max();
			found = false;
			contactPts.clear();
		}

		size_t numContactPoints() const {
			return contactPts.size();
		}
		FloatingType minDisToPrimitives;

	};

	// link face to VF contact stored by vertex
	struct FVContactInfo
	{
		IdType meshIdVSide = -1;
		IdType vertexId = -1;
		IdType contactId = -1;
	};

	struct EEContactPointInfo
	{
		int contactEdgeId1 = -1;
		int contactMeshId1 = -1;

		int contactEdgeId2 = -1;
		int contactMeshId2 = -1;

		FloatingType mu1 = 0.f;
		FloatingType mu2 = 0.f;

		// c1 = lerp(p1, p2, mu1); c2 = lerp(q1, q2, mu2); p1,p2 from edge 1 and q1,q2 from edge 2
		// lerp: 1.0f-t*v0 + t*v1
		Vec3 c1, c2; // c1 is from the current edge and c2 is from the opposite edge

		FloatingType d = 0.f;
	};

	struct ClothEEContactQueryResult {
		CPArray<EEContactPointInfo, EE_CONTACT_PREALLOCATE> contactPts;

		// query point info
		ClothContactDetector* pContactDetector = nullptr;

		int queryMeshId = -1;
		int queryPrimitiveId = -1;
		bool found = true;

		ClothEEContactQueryResult() {
			minDisToPrimitives = std::numeric_limits<FloatingType>::max();
		}

		void reset() {
			minDisToPrimitives = std::numeric_limits<FloatingType>::max();
			found = false;
			contactPts.clear();
		}

		size_t numContactPoints() const {
			return contactPts.size();
		}

		FloatingType minDisToPrimitives;
	};

	struct ClothContactDetector 
	{
		typedef std::shared_ptr<ClothContactDetector> SharedPtr;
		typedef ClothContactDetector* Ptr;

		ClothContactDetector(const ClothContactDetectorParameters::SharedPtr pParameters);
		void initialize(std::vector<TriMeshFEM::SharedPtr> in_targetMeshes);
		bool contactQueryVF(IdType meshId, IdType vId, ClothVFContactQueryResult* pResult);
		bool contactQueryFV(IdType meshId, IdType fId, ClothVFContactQueryResult* pResult, IdType centerVId=-1 );
		bool contactQueryEE(IdType meshId, IdType vId, ClothEEContactQueryResult* pResult);

		void updateBVH(RTCBuildQuality sceneQuality = RTC_BUILD_QUALITY_REFIT);

		void resetFaceContactInfo() 
		{
			for (size_t i = 0; i < faceMinDisToVertices.size(); i++){
				for (size_t j = 0; j < faceMinDisToVertices[i].size(); j++){
					faceMinDisToVertices[i][j] = pParams->maxQueryDis;
					faceContactInfos[i][j].clear();
				}
			}
		}

		CPArrayStaticAtomic<FVContactInfo, FV_CONTACT_PREALLOCATE>& getFaceContactInfo(IdType meshId, IdType faceId)
		{
			return faceContactInfos[meshId][faceId];
		}

		CFloatingType getFaceMinDis(IdType meshId, IdType faceId) { return faceMinDisToVertices[meshId][faceId]; }

		const ClothContactDetectorParameters& parameters() const
		{
			return *(ClothContactDetectorParameters*)pParams.get();
		}

	public:
		ClothContactDetectorParameters::SharedPtr pParams;

		RTCScene targetMeshEdgesScene;
		// only used when supportFVQuery is true
		RTCScene targetMeshVerticesScene;
		// numMesh x numFaces
		std::vector<std::vector<std::atomic<FloatingType>>> faceMinDisToVertices;
		std::vector<std::vector<CPArrayStaticAtomic<FVContactInfo, FV_CONTACT_PREALLOCATE>>> faceContactInfos;
		std::vector<TriMeshFEM::SharedPtr> targetMeshes;

		RTCScene targetMeshFacesScene;
		RTCDevice device;
	};

	void updateVFContactPointInfo(std::vector<std::shared_ptr<TriMeshFEM>>& meshPtrs, VFContactPointInfo& contactPointInfo);
	void updateEEContactPointInfo(std::vector<std::shared_ptr<TriMeshFEM>>& meshPtrs, EEContactPointInfo& contactPointInfo);
}
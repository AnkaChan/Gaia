#pragma once


#include <vector>
#include <memory>
#include <embree3/rtcore.h>

#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "../common/math/vec4.h"
#include "../common/math/affinespace.h"
#include "../common/math/constants.h"

#include "CollisionDetertionParameters.h"
#include "TriangleTriangleIntersection.h"

#define PREALLOCATED_NUM_TRI_TRI_COLLISIONS 16
#define PREALLOCATED_NUM_INCLUSIONS 8
#include "Types/Types.h"

#include <atomic>

namespace GAIA {
    struct TetMeshFEM;
	struct VolumetricCollisionDetector;

	struct TriTriIntersectionResults {
		void clear() {
			numCollisions = 0;
			if (overFlow) {
				maxNumTriTriIntersections = size_t(1.5 * intersections.size());
				intersections.resize(maxNumTriTriIntersections);
				overFlow = false;
			}
			collidedTriangles.clear();
		}

		size_t maxNumTriTriIntersections = -1;
		std::atomic<int> numCollisions = 0;
		bool overFlow = false;
		std::vector<GAIA::TriTriIntersection> intersections;

		std::vector<std::pair<IdType, IdType>> collidedTriangles;

	};

	struct PointInclusionRTCContext : public RTCIntersectContext
	{
		VolumetricCollisionDetector* pVolCol;
		int numHits = 0;
		// size is the number of meshes, each element record the number of ray intersections 
		VecDynamicI numRayIntersections;
		// whether the point intersects with and edge 
		bool sucess = true;
	};

	struct SurfaceMeshCloestPointQueryResult {
		Vec3 barys;
		Vec3 closestPt;
		Vec3I faceVIds;
		int intersectingMeshId = -1;
		int closestFaceId;
		VolumetricCollisionDetector* pCollisionDetector;
		ClosestPointOnTriangleType closestPointType;
	};


	struct TriangleCollisionResults {
		void clear()
		{
			numIntersections = 0;
			includedMeshes.clear();
			intersectionPoints.clear();
			includedVerts = 0;
		}
		
		std::atomic<int> numIntersections;
		// points to tritriIntersectionResults of VolumetricCollisionDetector
		int intersections[PREALLOCATED_NUM_TRI_TRI_COLLISIONS];
		// barycentrics
		CPArray<Vec3, PREALLOCATED_NUM_TRI_TRI_COLLISIONS> intersectionPoints;
		int includedVerts;

		CPArray<int, PREALLOCATED_NUM_INCLUSIONS> includedMeshes;

	};

	struct VolumetricCollisionDetector
	{
		VolumetricCollisionDetector(const VolumetricCollisionParameters& in_params);
		~VolumetricCollisionDetector();
		void initialize(std::vector<std::shared_ptr<TetMeshFEM>> tMeshes);
		void updateBVH( RTCBuildQuality sceneQuality);
		// Do not run this for multiple meshes in parallel, but it can be run in parallel for all the surface triangles* of a mesh
		
		// tri-tri intersection test
		void clearTriangleIntersections();
		void triangleIntersectionTest();
		void findTriangleIntersectionPolygons();
		void clusterIntersectionPoints(int meshId, int triId);

		// point inclusion test
		void clearpointInclusionResults();
		void closestSurfacePtQuery(int iMesh, int iV, int intersectingMeshId, SurfaceMeshCloestPointQueryResult* pResult);
		void pointInclusionTest();


		int numSurfaceFaces(int meshId);
		Vec3 getIntersectionPosition(int iIntersection, int triId, int meshId);


		std::vector<TriangleCollisionResults*> triangleIntersections;
		TriTriIntersectionResults tritriIntersectionResults;

		std::vector<std::vector<PointInclusionRTCContext>> pointInclusionTestResults;


		const VolumetricCollisionParameters& params;
		std::vector<std::shared_ptr<TetMeshFEM>> tMeshPtrs;
		//std::vector<std::vector<int>> pointsForInclusionTest;

		//std::vector<RTCScene> surfaceMeshScenes;
		RTCScene surfaceMeshScene = nullptr;
		RTCScene surfaceMeshSceneRayTracing = nullptr;
		RTCDevice device = nullptr;


	};


}
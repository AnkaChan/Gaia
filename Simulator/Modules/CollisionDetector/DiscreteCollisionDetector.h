#pragma once

// looks like this include file need to be included before others

#include <vector>
#include <memory>
#include <embree3/rtcore.h>

#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "../common/math/vec4.h"
#include "../common/math/affinespace.h"
#include "../common/math/constants.h"

#include "CollisionDetertionParameters.h"

#include "CollisionGeometry.h"

#include "../TetMesh/TetMeshFEM.h"

namespace GAIA {
    struct TetMeshFEM;
    struct DiscreteCollisionDetector;

    struct ClosestPointQueryResult
    {
        ClosestPointQueryResult()
        {}

        int idVQuery = -1;
        int idEmbraceTet = -1;
        int idTMQuery = -1;

        int closestFaceId = -1;
        embree::Vec3fa closestPt;
        embree::Vec3fa closestPtBarycentrics;
        ClosestPointOnTriangleType closestPointType;

        // if the closest point on the surface can be found
        bool found = false;

        bool checkFeasibleRegion = false;
        bool checkTetTraverse = true;

        DiscreteCollisionDetector* pDCD = nullptr;

        int numberOfBVHQuery = 0;
        int numberOfTetTraversal = 0;
        int numberOfTetsTraversed = 0;
        //int nFaceTraversed = 0;
        //PathFinder::TetMeshPathFinder* pathFinder = nullptr;
        //PathFinder::RayTargetPointIntersectionType intersectionType;
        //std::vector<PathFinder::TM::TPtr>* traversedTVec = nullptr;

        //void* pMeshClosestElement;

    };


	struct DiscreteCollisionDetector
	{
		DiscreteCollisionDetector(const CollisionDetectionParamters & in_params);
		void initialize(std::vector<std::shared_ptr<TetMeshFEM>> tMeshes);

        void updateBVH(RTCBuildQuality tetMeshSceneQuality, RTCBuildQuality surfaceSceneQuality
            , bool updateSurfaceScene);

        // vId: index of tetmesh vertex (not surface vertex, this also works for interior verts)
        bool vertexCollisionDetection(int32_t vId, int32_t tMeshId, CollisionDetectionResult* pResult);
        bool closestPointQuery(CollisionDetectionResult* pResult, ClosestPointQueryResult* pClosestPtResult);

        // edgeID: 0,1,2 represents 
        bool checkFeasibleRegion(embree::Vec3fa& p, TetMeshFEM *pTM, int32_t faceId, 
            ClosestPointOnTriangleType pointType, float feasibleREgionEpsilon);

		RTCScene tetMeshesScene;
		int numTetsTotal;

		std::vector<std::shared_ptr<TetMeshFEM>> tMeshPtrs;
		// a scene for each surface mesh
        // used for geodesic closest surface point query 
		std::vector<RTCScene> surfaceMeshScenes;
		
        void computeNormal(CollisionDetectionResult& colResult, int32_t iIntersection, Vec3& normal);

		RTCDevice device;

		const CollisionDetectionParamters& params;

	};

    embree::Vec3fa loadVertexPos(TetMeshFEM* pTM, int32_t vId);
    embree::Vec3fa faceNormal(TetMeshFEM* pTM, int32_t faceId);

    // edgeID: 0,1,2 represents AB, BC, CA respectively
    bool checkEdgeFeasibleRegion(embree::Vec3fa& p, TetMeshFEM* pTM, int32_t faceId, int32_t edgeId,
        int32_t edgeVId1, int32_t edgeVId2, float feasibleRegionEpsilon);

    // vId is a tet vertex id, but it must be a surface vertex, i.e., pTM->
    bool checkVertexFeasibleRegion(embree::Vec3fa& p, TetMeshFEM* pTM, int32_t vId, float feasibleRegionEpsilon);

    inline embree::Vec3fa GAIA::loadVertexPos(TetMeshFEM* pTM, int32_t vId)
    {
        return embree::Vec3fa::loadu(pTM->mVertPos.col(vId).data());
    }

    inline embree::Vec3fa GAIA::faceNormal(TetMeshFEM* pTM, int32_t faceId)
    {
        embree::Vec3fa a = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(0, faceId)),
            b = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(1, faceId)),
            c = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(2, faceId));

        embree::Vec3fa ab = b - a;
        embree::Vec3fa ac = c - a;

        embree::Vec3fa normal = embree::cross(ab, ac);
        normal = embree::normalize(normal);

        return normal;
    }

    inline embree::Vec3fa faceOrientedArea(TetMeshFEM* pTM, int32_t faceId)
    {
        embree::Vec3fa a = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(0, faceId)),
            b = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(1, faceId)),
            c = loadVertexPos(pTM, pTM->surfaceFacesTetMeshVIds()(2, faceId));

        embree::Vec3fa ab = b - a;
        embree::Vec3fa ac = c - a;

        embree::Vec3fa orientedArea = embree::cross(ab, ac);

        return orientedArea;
    }


}
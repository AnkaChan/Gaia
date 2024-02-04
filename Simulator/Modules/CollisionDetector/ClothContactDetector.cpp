#include "ClothContactDetector.h"
#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "../common/math/vec4.h"
#include "../common/math/affinespace.h"
#include "../common/math/constants.h"
#include "../CollisionDetector/CollisionGeometry.h"

#include "TriMeshCollisionGeometry.h"

//#define SKIP_FEASIBLE_REGION_CHECK

using namespace GAIA;

bool triMeshRadiusQueryWithTopologyFilteringFunc(RTCPointQueryFunctionArguments* args)
{
    ClothVFContactQueryResult* result = (ClothVFContactQueryResult*)args->userPtr;
    // TetMeshFEM* pTMQuery = result->pDCD->tMeshPtrs[result->idTMQuery].get();

    const TriMeshFEM* pMeshQuery = result->pMeshQuery;
    const int queryPremitiveId = result->queryPrimitiveId;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    const unsigned int primID = args->primID;
    MeshClosestPointQuery* pMatcher = result->pMatcher;
    TriMeshFEM* pTargetMesh = pMatcher->targetMeshes[geomID].get();

    if (pMeshQuery == pTargetMesh)
    {
        // filter out the adjacent face
        for (size_t faceNeiVId = 0; faceNeiVId < 3; faceNeiVId++)
        {
            int neiVId = pTargetMesh->facePosVId(primID, faceNeiVId);
            if (neiVId == queryPremitiveId)
            {
                return false;
            }
        }
    }


    const embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    //const embree::Vec3ia face(pTargetMesh->facePos(0, primID),
    //    pTargetMesh->facePos(1, primID), pTargetMesh->facePos(2, primID));
    const IdType* face = pTargetMesh->facePos.col(primID).data();

    const embree::Vec3fa a = embree::Vec3fa::loadu(pTargetMesh->vertex(face[0]).data());
    const embree::Vec3fa b = embree::Vec3fa::loadu(pTargetMesh->vertex(face[1]).data());
    const embree::Vec3fa c = embree::Vec3fa::loadu(pTargetMesh->vertex(face[2]).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    const embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(queryPt, closestP);

    if (d < args->query->radius)
    {
        // evalute whether this closest has been added 
        int primitiveId = -1;
        switch (pointType)
        {
        case GAIA::ClosestPointOnTriangleType::AtA:
            primitiveId = pTargetMesh->facePosVId(primID, 0);
            break;
        case GAIA::ClosestPointOnTriangleType::AtB:
            primitiveId = pTargetMesh->facePosVId(primID, 1);
            break;
        case GAIA::ClosestPointOnTriangleType::AtC:
            primitiveId = pTargetMesh->facePosVId(primID, 2);
            break;
        case GAIA::ClosestPointOnTriangleType::AtAB:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(0, primID);
            break;
        case GAIA::ClosestPointOnTriangleType::AtBC:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(1, primID);
            break;
        case GAIA::ClosestPointOnTriangleType::AtAC:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(2, primID);
            break;
        case GAIA::ClosestPointOnTriangleType::AtInterior:
            primitiveId = primID;
            break;
        case GAIA::ClosestPointOnTriangleType::NotFound:
            break;
        default:
            break;
        }
        
        ClosestPointOnPrimitiveType primitiveType = getClosestPointOnPrimitiveType(pointType);

        for (size_t iClosestP = 0; iClosestP < result->contactPts.size(); iClosestP++)
        {
            if (result->contactPts[iClosestP].primitiveType == primitiveType
                && result->contactPts[iClosestP].primitiveId == primitiveId)
            {
					return false;
			}
        }

#ifndef SKIP_FEASIBLE_REGION_CHECK

        bool inFeasibleRegion = checkFeasibleRegion(queryPt, pTargetMesh, primID, pointType, 1e-3);
        if (!inFeasibleRegion)
        {
            return false;
        }
#endif // !SKIP_FEASIBLE_REGION_CHECK
        
        result->contactPts.emplace_back();

        result->contactPts.back().closestFaceId = primID;
        result->contactPts.back().closestMeshId = geomID;
        result->contactPts.back().d = d;
        result->contactPts.back().closestPt << closestP.x, closestP.y, closestP.z;
        result->contactPts.back().barycentrics << closestPtBarycentrics.x, closestPtBarycentrics.y, closestPtBarycentrics.z;
        result->contactPts.back().closestPtType = pointType;
       
        result->contactPts.back().primitiveId = primitiveId;
        result->contactPts.back().primitiveType = primitiveType;

        //result->closestPtBarycentrics = closestPtBarycentrics;

        // record that at least one closest point search has succeeded
        result->found = true;

        return false; // Return true to indicate that the query radius changed.
    }

    return false;
}

bool triMeshEdgeContactQueryFunc(RTCPointQueryFunctionArguments* args)
{
    ClothEEContactQueryResult* result = (ClothEEContactQueryResult*)args->userPtr;
    const TriMeshFEM* pMeshQuery = result->pMeshQuery;
    const int queryPremitiveId = result->queryPrimitiveId;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    const unsigned int primID = args->primID;
    ClothContactDetector* pContactDetector = result->pContactDetector;
    TriMeshFEM* pTargetMesh = pContactDetector->targetMeshes[geomID].get();

    const EdgeInfo& edgeInfoQuery = pMeshQuery->pTopology->edgeInfos[queryPremitiveId];
    const EdgeInfo& edgeInfoTarget = pTargetMesh->pTopology->edgeInfos[primID];

    if (pMeshQuery == pTargetMesh)
    {
        if (primID == queryPremitiveId)
        {
            return false;
        }
        if (edgeInfoQuery.eV1 == edgeInfoTarget.eV1 
            || edgeInfoQuery.eV1 == edgeInfoTarget.eV2
            || edgeInfoQuery.eV2 == edgeInfoTarget.eV1
			|| edgeInfoQuery.eV2 == edgeInfoTarget.eV2)
        {
            return false;
        }
    }


    const embree::Vec3fa p1 = embree::Vec3fa::loadu(pMeshQuery->vertex(edgeInfoQuery.eV1).data());
    const embree::Vec3fa p2 = embree::Vec3fa::loadu(pMeshQuery->vertex(edgeInfoQuery.eV2).data());

    const embree::Vec3fa q1 = embree::Vec3fa::loadu(pMeshQuery->vertex(edgeInfoTarget.eV1).data());
    const embree::Vec3fa q2 = embree::Vec3fa::loadu(pMeshQuery->vertex(edgeInfoTarget.eV2).data());
    embree::Vec3fa c1, c2;
    FloatingType mua, mub;
    get_closest_points_between_segments(p1, p2, q1, q2, c1, c2, mua, mub);

    FloatingType d = embree::distance(c1, c2);

    if (d < pContactDetector->pParams->maxQueryDis
        && mua > 0.f && mua < 1.f
        && mub > 0.f && mub < 1.f
        ) // otherwise it degenerates to a v-f contact case
    {
        result->found = true;

        result->contactPts.emplace_back();
        result->contactPts.back().closestEdgeId = primID;
        result->contactPts.back().closestMeshId = geomID;
        result->contactPts.back().mu_this = mua;
        result->contactPts.back().mu_opposite = mub;
        result->contactPts.back().c1 << c1.x, c1.y, c1.z;
        result->contactPts.back().c2 << c2.x, c2.y, c2.z;
        result->contactPts.back().d = d;
    }

    return false;
}

void edgeBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const TriMeshFEM* pMesh = (const TriMeshFEM*)args->geometryUserPtr;

    embree::BBox3fa bounds = embree::empty;

    const EdgeInfo& edgeInfo = pMesh->pTopology->edgeInfos[args->primID];
    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->vertex(edgeInfo.eV1).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->vertex(edgeInfo.eV2).data());
    bounds.extend(a);

    *(embree::BBox3fa*)args->bounds_o = bounds;
}

GAIA::ClothContactDetector::ClothContactDetector(const ClothContactDetectorParameters::SharedPtr in_pParams)
	: MeshClosestPointQuery(in_pParams)
{

}

void GAIA::ClothContactDetector::initialize(std::vector<TriMeshFEM::SharedPtr> in_targetMeshes)
{
	MeshClosestPointQuery::initialize(in_targetMeshes);

    for (size_t meshId = 0; meshId < targetMeshes.size(); meshId++)
    {
        RTCGeometry geom = rtcGetGeometry(targetMeshScene, meshId);

        rtcSetGeometryPointQueryFunction(geom, triMeshRadiusQueryWithTopologyFilteringFunc);
        rtcCommitGeometry(geom);
    }
    rtcSetSceneFlags(targetMeshScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(targetMeshScene, RTC_BUILD_QUALITY_LOW);
    rtcCommitScene(targetMeshScene);

    targetMeshEdgesScene = rtcNewScene(device);
    rtcSetSceneFlags(targetMeshEdgesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(targetMeshScene, RTC_BUILD_QUALITY_LOW);

    for (size_t meshId = 0; meshId < targetMeshes.size(); meshId++)
    {
		TriMeshFEM* pMesh = targetMeshes[meshId].get();

        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geom, pMesh->numEdges());
        rtcSetGeometryBoundsFunction(geom, edgeBoundsFunc, nullptr);

        rtcSetGeometryUserData(geom, (void*)(pMesh));

        rtcSetGeometryPointQueryFunction(geom, triMeshEdgeContactQueryFunc);
        rtcCommitGeometry(geom);
        unsigned int geomId = meshId;
        rtcAttachGeometryByID(targetMeshEdgesScene, geom, geomId);
        rtcReleaseGeometry(geom);
	}


}

bool GAIA::ClothContactDetector::contactQueryVF(TriMeshFEM* pClothMesh, IdType vId, ClothVFContactQueryResult* pResult)
{
    pResult->reset();
    RTCPointQuery query;

	Vec3 p = pClothMesh->vertex(vId);

    pResult->pMeshQuery = pClothMesh;
    pResult->queryPrimitiveId = vId;
    pResult->pMatcher = this;
        
    query.x = p(0);
    query.y = p(1);
    query.z = p(2);

    query.radius = pParams->maxQueryDis;
    query.time = 0.f;

    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(targetMeshScene, &query, &context, nullptr, (void*)pResult);

    return pResult->found;
}

bool GAIA::ClothContactDetector::contactQueryEE(TriMeshFEM* pClothMesh, IdType eId, ClothEEContactQueryResult* pResult)
{
    pResult->reset();
    RTCPointQuery query;

    const EdgeInfo& eInfo = pClothMesh->pTopology->edgeInfos[eId];

    FloatingType eLen = (pClothMesh->vertex(eInfo.eV1) - pClothMesh->vertex(eInfo.eV2)).norm();

    // query center will be the edge center
    Vec3 p = (pClothMesh->vertex(eInfo.eV1) + pClothMesh->vertex(eInfo.eV2)) * 0.5f;
    // query radius is the edge length + contact radius
    query.radius = (pParams->maxQueryDis + eLen) * 0.5f;
    query.time = 0.f;

    pResult->pMeshQuery = pClothMesh;
    pResult->queryPrimitiveId = eId;
    pResult->pContactDetector = this;

    query.x = p(0);
    query.y = p(1);
    query.z = p(2);
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(targetMeshEdgesScene, &query, &context, nullptr, (void*)pResult);

    return false;
}

void GAIA::ClothContactDetector::updateBVH(RTCBuildQuality sceneQuality)
{
    MeshClosestPointQuery::updateBVH(sceneQuality);

    RTCBuildQuality geomQuality = sceneQuality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }

    rtcSetSceneBuildQuality(targetMeshEdgesScene, sceneQuality);

for (size_t meshId = 0; meshId < targetMeshes.size(); meshId++)
	{
        RTCGeometry geom = rtcGetGeometry(targetMeshEdgesScene, meshId);
        rtcSetGeometryBuildQuality(geom, geomQuality);
        rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geom);
	}

	rtcCommitScene(targetMeshEdgesScene);
}

bool GAIA::ClothContactDetectorParameters::fromJson(nlohmann::json& j)
{
    MeshClosestPointQueryParameters::fromJson(j);

    return true;
}

bool GAIA::ClothContactDetectorParameters::toJson(nlohmann::json& j)
{
    MeshClosestPointQueryParameters::toJson(j);
    return true;
}

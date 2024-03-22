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

// if the closest point is on the edge, the contact will be detected twice, so is the case when the closest point is on the vertex,
// we need to avoid those duplications; this is done by recording the id of the contact primitive and the type of the closest point
// when a new closest point is found, we check whether it has been added before
bool triMeshVFRadiusQueryWithTopologyFilteringFunc(RTCPointQueryFunctionArguments* args)
{
    ClothVFContactQueryResult* result = (ClothVFContactQueryResult*)args->userPtr;
    // TetMeshFEM* pTMQuery = result->pDCD->tMeshPtrs[result->idTMQuery].get();
    assert(args->userPtr);

    ClothContactDetector* pContactDetector = result->pContactDetector;
    IdType queryMeshId = result->queryMeshId;
    const TriMeshFEM* pMeshQuery = pContactDetector->targetMeshes[queryMeshId].get();
    // face Id
    const int queryPremitiveId_vertex = result->queryPrimitiveId;

    const unsigned int geomID_face = args->geomID;
    const unsigned int primID_face = args->primID;
    TriMeshFEM* pTargetMesh = pContactDetector->targetMeshes[geomID_face].get();

    if (pMeshQuery == pTargetMesh)
    {
        // filter out the adjacent face
        for (size_t faceNeiVId = 0; faceNeiVId < 3; faceNeiVId++)
        {
            int neiVId = pTargetMesh->facePosVId(primID_face, faceNeiVId);
            if (neiVId == queryPremitiveId_vertex)
            {
                return false;
            }
        }
    }

    const embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    //const embree::Vec3ia face(pTargetMesh->facePos(0, primID),
    //    pTargetMesh->facePos(1, primID), pTargetMesh->facePos(2, primID));
    const IdType* face = pTargetMesh->facePos.col(primID_face).data();

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
            primitiveId = pTargetMesh->facePosVId(primID_face, 0);
            break;
        case GAIA::ClosestPointOnTriangleType::AtB:
            primitiveId = pTargetMesh->facePosVId(primID_face, 1);
            break;
        case GAIA::ClosestPointOnTriangleType::AtC:
            primitiveId = pTargetMesh->facePosVId(primID_face, 2);
            break;
        case GAIA::ClosestPointOnTriangleType::AtAB:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(0, primID_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtBC:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(1, primID_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtAC:
            primitiveId = pTargetMesh->pTopology->faces3NeighborEdges(2, primID_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtInterior:
            primitiveId = primID_face;
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

        bool inFeasibleRegion = checkFeasibleRegion(queryPt, pTargetMesh, primID_face, pointType, 1e-3);
        if (!inFeasibleRegion)
        {
            return false;
        }
#endif // !SKIP_FEASIBLE_REGION_CHECK
        
        result->contactPts.emplace_back();

        result->contactPts.back().contactVertexId = queryPremitiveId_vertex;
        result->contactPts.back().contactVertexSideMeshId = queryMeshId;
        result->contactPts.back().contactFaceId = primID_face;
        result->contactPts.back().contactFaceSideMeshId = geomID_face;
        result->contactPts.back().d = d;
        result->contactPts.back().contactPoint << closestP.x, closestP.y, closestP.z;

        computeVFContactNormalTriMesh(a, b, c, queryPt, closestP, pointType, result->contactPts.back().contactPointNormal);

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

bool triMeshFVRadiusQueryWithTopologyFilteringFunc(RTCPointQueryFunctionArguments* args)
{
    ClothVFContactQueryResult* result = (ClothVFContactQueryResult*)args->userPtr;

    ClothContactDetector* pContactDetector = result->pContactDetector;
    IdType queryMeshId_face = result->queryMeshId;
    const TriMeshFEM* pMeshQuery_faceSide = pContactDetector->targetMeshes[queryMeshId_face].get();
    const int queryPremitiveId_face = result->queryPrimitiveId;

    assert(args->userPtr);
    const unsigned int geomID_vertex = args->geomID;
    const unsigned int primID_vertex = args->primID;
    IdType vId = primID_vertex;
    TriMeshFEM* pTargetMesh = pContactDetector->targetMeshes[geomID_vertex].get();

    const IdType* face = pMeshQuery_faceSide->facePos.col(queryPremitiveId_face).data();

    if (geomID_vertex == queryMeshId_face)
    {
        // filter out the vertices from the same face
        for (size_t iFV = 0; iFV < 3; iFV++)
        {
            int faceVId = face[iFV];
            if (faceVId == vId)
            {
                return false;
            }
        }
    }

    const embree::Vec3fa v = embree::Vec3fa::loadu(pTargetMesh->vertex(vId).data());

    //const embree::Vec3ia face(pTargetMesh->facePos(0, primID),
    //    pTargetMesh->facePos(1, primID), pTargetMesh->facePos(2, primID));

    const embree::Vec3fa a = embree::Vec3fa::loadu(pMeshQuery_faceSide->vertex(face[0]).data());
    const embree::Vec3fa b = embree::Vec3fa::loadu(pMeshQuery_faceSide->vertex(face[1]).data());
    const embree::Vec3fa c = embree::Vec3fa::loadu(pMeshQuery_faceSide->vertex(face[2]).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    const embree::Vec3fa closestP = GAIA::closestPointTriangle(v, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(v, closestP);
    
    float contactRadius = pContactDetector->parameters().maxQueryDis;

    if (d < contactRadius)
    {
        // evalute whether this closest has been added 
        int primitiveId = -1;
        switch (pointType)
        {
        // if the closest point is on the vertex, it is either the center vertex, which will be handled by the VF query
        // or it's one of the vertices of the face, which will result in zero force on the center vertex
        case GAIA::ClosestPointOnTriangleType::AtA:
            // primitiveId = pMeshQuery_faceSide->facePosVId(queryPremitiveId_face, 0);
        case GAIA::ClosestPointOnTriangleType::AtB:
            // primitiveId = pMeshQuery_faceSide->facePosVId(queryPremitiveId_face, 1);
        case GAIA::ClosestPointOnTriangleType::AtC:
            // primitiveId = pMeshQuery_faceSide->facePosVId(queryPremitiveId_face, 2);
            return false;
            break;
        case GAIA::ClosestPointOnTriangleType::AtAB:
            primitiveId = pMeshQuery_faceSide->pTopology->faces3NeighborEdges(0, queryPremitiveId_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtBC:
            primitiveId = pMeshQuery_faceSide->pTopology->faces3NeighborEdges(1, queryPremitiveId_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtAC:
            primitiveId = pMeshQuery_faceSide->pTopology->faces3NeighborEdges(2, queryPremitiveId_face);
            break;
        case GAIA::ClosestPointOnTriangleType::AtInterior:
            primitiveId = queryPremitiveId_face;
            break;
        case GAIA::ClosestPointOnTriangleType::NotFound:
            break;
        default:
            break;
        }

        ClosestPointOnPrimitiveType primitiveType = getClosestPointOnPrimitiveType(pointType);
        // skip the collision for the center vertex since it's already handled by the vertex-face contact query
        //if (primitiveType == ClosestPointOnPrimitiveType::AtVertex && result->centerVertexId == primitiveId)
        //{
        //    return false;
        //}
        //for (size_t iClosestP = 0; iClosestP < result->contactPts.size(); iClosestP++)
        //{
        //    if (result->contactPts[iClosestP].primitiveType == primitiveType
        //        && result->contactPts[iClosestP].primitiveId == primitiveId)
        //    {
        //        return false;
        //    }
        //}

#ifndef SKIP_FEASIBLE_REGION_CHECK

        bool inFeasibleRegion = checkFeasibleRegion(v, pMeshQuery_faceSide, queryPremitiveId_face, pointType, 1e-3);
        if (!inFeasibleRegion)
        {
            return false;
        }
#endif // !SKIP_FEASIBLE_REGION_CHECK

        result->contactPts.emplace_back();

        // contary to the VF query, the contact face id is queryPremitiveId
        result->contactPts.back().contactVertexId = primID_vertex;
        result->contactPts.back().contactVertexSideMeshId = geomID_vertex;
        result->contactPts.back().contactFaceId = queryPremitiveId_face;
        result->contactPts.back().contactFaceSideMeshId = queryMeshId_face;

        result->contactPts.back().contactPoint << closestP.x, closestP.y, closestP.z;

        computeVFContactNormalTriMesh(a, b, c, v, closestP, pointType, result->contactPts.back().contactPointNormal);

        result->contactPts.back().d = d;
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

// this function is used to detect edge-edge contacts
// the contact can only be detected when the closest point is on the interior of the edge
// if the closest point is on the vertex, it's already handled by the vertex-face contact query
// therefore, no duplicate contact will be detected
bool triMeshEdgeContactQueryFunc(RTCPointQueryFunctionArguments* args)
{
    ClothEEContactQueryResult* result = (ClothEEContactQueryResult*)args->userPtr;
    ClothContactDetector* pContactDetector = result->pContactDetector;
    IdType queryMeshId = result->queryMeshId;
    const TriMeshFEM* pMeshQuery = pContactDetector->targetMeshes[queryMeshId].get();
    const int queryPremitiveId = result->queryPrimitiveId;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    const unsigned int primID = args->primID;
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
        result->contactPts.back().contactEdgeId1 = queryPremitiveId;
        result->contactPts.back().contactMeshId1 = queryMeshId;

        result->contactPts.back().contactEdgeId2 = primID;
        result->contactPts.back().contactMeshId2 = geomID;
        result->contactPts.back().mu1 = mua;
        result->contactPts.back().mu2 = mub;
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

void verticesBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const TriMeshFEM* pMesh = (const TriMeshFEM*)args->geometryUserPtr;

    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->vertex(args->primID).data());
    embree::BBox3fa& bounds = *(embree::BBox3fa*)args->bounds_o;
    bounds.upper = bounds.lower = a;
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
        RTCGeometry geom = rtcGetGeometry(targetMeshFacesScene, meshId);

        rtcSetGeometryPointQueryFunction(geom, triMeshVFRadiusQueryWithTopologyFilteringFunc);
        rtcCommitGeometry(geom);
    }

    // change the call back function for the targetMeshFacesScene to the VF query contact function
    rtcSetSceneFlags(targetMeshFacesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(targetMeshFacesScene, RTC_BUILD_QUALITY_LOW);
    rtcCommitScene(targetMeshFacesScene);

    targetMeshEdgesScene = rtcNewScene(device);
    rtcSetSceneFlags(targetMeshEdgesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(targetMeshEdgesScene, RTC_BUILD_QUALITY_LOW);

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
    rtcCommitScene(targetMeshEdgesScene);


    if (parameters().supportFVQuery)
    {
        targetMeshVerticesScene = rtcNewScene(device);
        rtcSetSceneFlags(targetMeshVerticesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
        rtcSetSceneBuildQuality(targetMeshVerticesScene, RTC_BUILD_QUALITY_LOW);

        for (size_t meshId = 0; meshId < targetMeshes.size(); meshId++)
        {
            TriMeshFEM* pMesh = targetMeshes[meshId].get();

            RTCGeometry geomVerts = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
            rtcSetGeometryUserPrimitiveCount(geomVerts, pMesh->numVertices());
            rtcSetGeometryBoundsFunction(geomVerts, verticesBoundsFunc, nullptr);

            rtcSetGeometryUserData(geomVerts, (void*)(pMesh));

            rtcSetGeometryPointQueryFunction(geomVerts, triMeshFVRadiusQueryWithTopologyFilteringFunc);
            rtcCommitGeometry(geomVerts);
            unsigned int geomId = meshId;
            rtcAttachGeometryByID(targetMeshVerticesScene, geomVerts, geomId);
            rtcReleaseGeometry(geomVerts);
        }
        rtcCommitScene(targetMeshVerticesScene);
    }

}

bool GAIA::ClothContactDetector::contactQueryVF(IdType meshId, IdType vId, ClothVFContactQueryResult* pResult)
{
    pResult->reset();
    RTCPointQuery query;
    TriMeshFEM* pClothMesh = targetMeshes[meshId].get();
	Vec3 p = pClothMesh->vertex(vId);

    pResult->queryMeshId = meshId;
    pResult->queryPrimitiveId = vId;
    pResult->pContactDetector = this;
        
    query.x = p(0);
    query.y = p(1);
    query.z = p(2);

    query.radius = pParams->maxQueryDis;
    query.time = 0.f;

    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(targetMeshFacesScene, &query, &context, nullptr, (void*)pResult);

    return pResult->found;
}

bool GAIA::ClothContactDetector::contactQueryFV(IdType meshId, IdType fId, ClothVFContactQueryResult* pResult, IdType centerVId)
{
    pResult->reset();
    RTCPointQuery query;
    TriMeshFEM* pClothMesh = targetMeshes[meshId].get();

    IdType * face = pClothMesh->facePos.col(fId).data();
    // query center will be the face's center
    Vec3 p = (pClothMesh->vertex(face[0]) + pClothMesh->vertex(face[1]) + pClothMesh->vertex(face[2])) / 3.f;

    // query radius is the radius of the triangle's circumcircle + contact radius
    query.radius = (p - pClothMesh->vertex(face[0])).norm() + pParams->maxQueryDis;
    query.time = 0.f;

    pResult->queryMeshId = meshId;
    pResult->queryPrimitiveId = fId;
    pResult->pContactDetector = this;
    // pResult->centerVertexId = centerVId;

    query.x = p(0);
    query.y = p(1);
    query.z = p(2);
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(targetMeshVerticesScene, &query, &context, nullptr, (void*)pResult);
    return false;
}

bool GAIA::ClothContactDetector::contactQueryEE(IdType meshId, IdType eId, ClothEEContactQueryResult* pResult)
{
    pResult->reset();
    RTCPointQuery query;
    TriMeshFEM* pClothMesh = targetMeshes[meshId].get();

    const EdgeInfo& eInfo = pClothMesh->pTopology->edgeInfos[eId];

    FloatingType eLen = (pClothMesh->vertex(eInfo.eV1) - pClothMesh->vertex(eInfo.eV2)).norm();

    // query center will be the edge center
    Vec3 p = (pClothMesh->vertex(eInfo.eV1) + pClothMesh->vertex(eInfo.eV2)) * 0.5f;
    // query radius is the edge length + contact radius
    query.radius = pParams->maxQueryDis + eLen * 0.5f;
    query.time = 0.f;

    pResult->queryMeshId = meshId;
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

    if (parameters().supportFVQuery)
    {
        rtcSetSceneBuildQuality(targetMeshVerticesScene, sceneQuality);
        for (size_t meshId = 0; meshId < targetMeshes.size(); meshId++)
        {
            RTCGeometry geom = rtcGetGeometry(targetMeshVerticesScene, meshId);
            rtcSetGeometryBuildQuality(geom, geomQuality);
            rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
            rtcCommitGeometry(geom);
        }
        rtcCommitScene(targetMeshVerticesScene);
    }

}

bool GAIA::ClothContactDetectorParameters::fromJson(nlohmann::json& j)
{
    MeshClosestPointQueryParameters::fromJson(j);
    EXTRACT_FROM_JSON(j, supportFVQuery);
    return true;
}

bool GAIA::ClothContactDetectorParameters::toJson(nlohmann::json& j)
{
    MeshClosestPointQueryParameters::toJson(j);
    PUT_TO_JSON(j, supportFVQuery);
    return true;
}

void GAIA::updateVFContactPointInfo(std::vector<std::shared_ptr<TriMeshFEM>>& meshPtrs, VFContactPointInfo& contactPointInfo)
{
    IdType meshIdVertexSide = contactPointInfo.contactVertexSideMeshId;
    IdType vertexId = contactPointInfo.contactVertexSideMeshId; 
    IdType meshIdFaceSide = contactPointInfo.contactFaceSideMeshId;
    IdType faceId = contactPointInfo.contactVertexId;

    const TriMeshFEM* pMeshVertexSide = meshPtrs[meshIdVertexSide].get();
    const TriMeshFEM* pMeshFaceSide = meshPtrs[meshIdFaceSide].get();

    const IdType* face = pMeshFaceSide->facePos.col(faceId).data();

    const embree::Vec3fa a = embree::Vec3fa::loadu(pMeshFaceSide->vertex(face[0]).data());
    const embree::Vec3fa b = embree::Vec3fa::loadu(pMeshFaceSide->vertex(face[1]).data());
    const embree::Vec3fa c = embree::Vec3fa::loadu(pMeshFaceSide->vertex(face[2]).data());

    const embree::Vec3fa p = embree::Vec3fa::loadu(pMeshVertexSide->vertex(vertexId).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    const embree::Vec3fa closestP = GAIA::closestPointTriangle(p, a, b, c, closestPtBarycentrics,
        pointType);
    ClosestPointOnPrimitiveType primitiveType = getClosestPointOnPrimitiveType(pointType);

    contactPointInfo.closestPtType = pointType;
    contactPointInfo.d = embree::distance(p, closestP);
    contactPointInfo.contactPoint << closestP.x, closestP.y, closestP.z;
    contactPointInfo.barycentrics << closestPtBarycentrics.x, closestPtBarycentrics.y, closestPtBarycentrics.z;

    computeVFContactNormalTriMesh(a, b, c, p, closestP, pointType, contactPointInfo.contactPointNormal);

};

void  GAIA::updateEEContactPointInfo(std::vector<std::shared_ptr<TriMeshFEM>>& meshPtrs, EEContactPointInfo& contactPointInfo)
{
    const TriMeshFEM* pMesh1 = meshPtrs[contactPointInfo.contactMeshId1].get();
    const TriMeshFEM* pMesh2 = meshPtrs[contactPointInfo.contactMeshId2].get();

    const EdgeInfo& edgeInfo1 = pMesh1->pTopology->edgeInfos[contactPointInfo.contactEdgeId1];
    const EdgeInfo& edgeInfo2 = pMesh2->pTopology->edgeInfos[contactPointInfo.contactEdgeId2];

    const embree::Vec3fa p1 = embree::Vec3fa::loadu(pMesh1->vertex(edgeInfo1.eV1).data());
    const embree::Vec3fa p2 = embree::Vec3fa::loadu(pMesh1->vertex(edgeInfo1.eV2).data());

    const embree::Vec3fa q1 = embree::Vec3fa::loadu(pMesh2->vertex(edgeInfo2.eV1).data());
    const embree::Vec3fa q2 = embree::Vec3fa::loadu(pMesh2->vertex(edgeInfo2.eV2).data());
    embree::Vec3fa c1, c2;
    FloatingType mua, mub;
    get_closest_points_between_segments(p1, p2, q1, q2, c1, c2, mua, mub);

    contactPointInfo.mu1 = mua;
    contactPointInfo.mu2 = mub;
    contactPointInfo.c1 << c1.x, c1.y, c1.z;
    contactPointInfo.c2 << c2.x, c2.y, c2.z;
}

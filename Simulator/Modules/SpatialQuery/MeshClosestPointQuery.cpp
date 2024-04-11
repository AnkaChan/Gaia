#include "MeshClosestPointQuery.h"
#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "../common/math/vec4.h"
#include "../common/math/affinespace.h"
#include "../common/math/constants.h"
#include "../CollisionDetector/CollisionGeometry.h"

using namespace GAIA;

bool triMeshClosestPointQueryFunc(RTCPointQueryFunctionArguments* args)
{
    TriMeshClosestPointQueryResult* result = (TriMeshClosestPointQueryResult*)args->userPtr;
    // TetMeshFEM* pTMQuery = result->pDCD->tMeshPtrs[result->idTMQuery].get();

    MeshClosestPointQuery* pMatcher = result->pMatcher;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    const unsigned int primID = args->primID;
    TriMeshFEM* pTargetMesh = pMatcher->targetMeshes[geomID].get();

    const embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    const embree::Vec3ia face(pTargetMesh->facePos(0, primID),
        pTargetMesh->facePos(1, primID), pTargetMesh->facePos(2, primID));

    const embree::Vec3fa a = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[0]).data());
    const embree::Vec3fa b = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[1]).data());
    const embree::Vec3fa c = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[2]).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    const embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(queryPt, closestP);

    if (d < args->query->radius)
    {
        if (result->contactPts.size() ==0)
        {
            result->contactPts.emplace_back();
            // record that at least one closest point search has succeeded
            result->found = true;
        }
        args->query->radius = d;

        result->contactPts.back().closestFaceId = primID;
        result->contactPts.back().d = d;
        result->contactPts.back().closestPt << closestP.x, closestP.y, closestP.z;
        result->contactPts.back().closestPtType = pointType;
        if (result->computeNormal)
        {
            result->contactPts.back().closestPtNormal = pTargetMesh->computeNormal(primID);
        }
        //result->closestPtBarycentrics = closestPtBarycentrics;

        return true; // Return true to indicate that the query radius changed.
    }

    return false;
}

bool triMeshRadiusQueryFunc(RTCPointQueryFunctionArguments* args)
{
    TriMeshClosestPointQueryResult* result = (TriMeshClosestPointQueryResult*)args->userPtr;
    // TetMeshFEM* pTMQuery = result->pDCD->tMeshPtrs[result->idTMQuery].get();

    MeshClosestPointQuery* pMatcher = result->pMatcher;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    const unsigned int primID = args->primID;
    TriMeshFEM* pTargetMesh = pMatcher->targetMeshes[geomID].get();

    const embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    const embree::Vec3ia face(pTargetMesh->facePos(0, primID),
        pTargetMesh->facePos(1, primID), pTargetMesh->facePos(2, primID));

    const embree::Vec3fa a = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[0]).data());
    const embree::Vec3fa b = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[1]).data());
    const embree::Vec3fa c = embree::Vec3fa::loadu(pTargetMesh->positions().col(face[2]).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    const embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(queryPt, closestP);

    if (d < args->query->radius)
    {
        result->contactPts.emplace_back();

        result->contactPts.back().closestFaceId = primID;
        result->contactPts.back().d = d;
        result->contactPts.back().closestPt << closestP.x, closestP.y, closestP.z;
        result->contactPts.back().closestPtType = pointType;

        if (result->computeNormal)
        {
            result->contactPts.back().closestPtNormal = pTargetMesh->computeNormal(primID);
        }

        //result->closestPtBarycentrics = closestPtBarycentrics;

        // record that at least one closest point search has succeeded
        result->found = true;

        return true; // Return true to indicate that the query radius changed.
    }

    return false;
}


GAIA::MeshClosestPointQuery::MeshClosestPointQuery(const MeshClosestPointQueryParameters::SharedPtr in_pParams)
	: pParams(in_pParams)
{
}

void GAIA::MeshClosestPointQuery::updateBVH(RTCBuildQuality sceneQuality)
{
    RTCBuildQuality geomQuality = sceneQuality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }

    rtcSetSceneBuildQuality(targetMeshFacesScene, sceneQuality);

    for (size_t iMesh = 0; iMesh < targetMeshes.size(); iMesh++)
    {
        // get the tet geom buffer
        if (targetMeshes[iMesh]->updated)
        {
            unsigned int geoId = iMesh;
            RTCGeometry geom = rtcGetGeometry(targetMeshFacesScene, geoId);

            rtcSetGeometryBuildQuality(geom, geomQuality);
            //rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
            //rtcCommitGeometry(geom);

            rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
            rtcCommitGeometry(geom);
        }
    }

    rtcCommitScene(targetMeshFacesScene);
}

void GAIA::MeshClosestPointQuery::initialize(std::vector<TriMeshFEM::SharedPtr> in_pTargetMesh)
{
	targetMeshes = in_pTargetMesh;

    device = rtcNewDevice(NULL);

    targetMeshFacesScene = rtcNewScene(device);
    rtcSetSceneFlags(targetMeshFacesScene, RTC_SCENE_FLAG_ROBUST);

    for (size_t iMesh = 0; iMesh < targetMeshes.size(); iMesh++)
    {
        RTCGeometry geomRTC = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

        rtcSetSceneBuildQuality(targetMeshFacesScene, RTC_BUILD_QUALITY_MEDIUM);

        rtcSetSharedGeometryBuffer(geomRTC,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, targetMeshes[iMesh]->positions().data(), 0, 3 * sizeof(float),
            targetMeshes[iMesh]->numVertices());

        rtcSetSharedGeometryBuffer(geomRTC,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, targetMeshes[iMesh]->facePos.data(), 0, 3 * sizeof(unsigned),
            targetMeshes[iMesh]->numFaces());

        if (pParams->onlyKeepClosest)
        {
            rtcSetGeometryPointQueryFunction(geomRTC, triMeshClosestPointQueryFunc);
        }
        else
        {
            rtcSetGeometryPointQueryFunction(geomRTC, triMeshRadiusQueryFunc);
        }

        rtcCommitGeometry(geomRTC);
        rtcAttachGeometryByID(targetMeshFacesScene, geomRTC, iMesh);
        rtcReleaseGeometry(geomRTC);
    }

    rtcCommitScene(targetMeshFacesScene);
}

bool GAIA::MeshClosestPointQuery::closestPointQuery(const Vec3 p, TriMeshClosestPointQueryResult* pClosestPtResult, bool computeNormal)
{
    pClosestPtResult->reset();
    pClosestPtResult->computeNormal = computeNormal;
    pClosestPtResult->pMatcher = this;
    RTCPointQuery query;

    query.x = p(0);
    query.y = p(1);
    query.z = p(2);

    query.radius = pParams->maxQueryDis;
    query.time = 0.f;

    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(targetMeshFacesScene, &query, &context, nullptr, (void*)pClosestPtResult);

    return pClosestPtResult->found;
}

bool GAIA::MeshClosestPointQueryParameters::fromJson(nlohmann::json& j)
{
    EXTRACT_FROM_JSON(j, maxQueryDis);
    return true;
}

bool GAIA::MeshClosestPointQueryParameters::toJson(nlohmann::json& j)
{
    PUT_TO_JSON(j, maxQueryDis);
    return true;
}

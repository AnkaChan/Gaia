#include "VolumetricCollisionDetector.h"
#include "CuMatrix/MatrixOps/CuMatrix.h"
#include "CuMatrix/Geometry/Geometry.h"
#include "../TetMesh/TetMeshFEM.h"
#include "CollisionGeometry.h"

#include <math.h>
#include "TriangleTriangleIntersection.h"
#include "../Parallelization/CPUParallelization.h"

#include <embree3/rtcore_common.h>

#define MAX_TOTAL_HITS 256
#define RAY_TRACING_NEAR 1e-4
#define RAY_EDGE_INTERSECTION_EPISLON 1e-4

using namespace GAIA;

void countAllHits(const struct RTCFilterFunctionNArguments* args);

void triangle_bounds_func(const struct RTCBoundsFunctionArguments* args)
{
    GAIA::TetMeshFEM* pMesh = (GAIA::TetMeshFEM*)args->geometryUserPtr;
    embree::BBox3fa bounds = embree::empty;

    const IdType* fVIds = pMesh->surfaceFacesTetMeshVIds().col(args->primID).data();
    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->mVertPos.col(fVIds[0]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->mVertPos.col(fVIds[1]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->mVertPos.col(fVIds[2]).data());
    bounds.extend(a);
   
    *(embree::BBox3fa*)args->bounds_o = bounds;
}

bool surfaceMeshClosestPointQueryFunc(RTCPointQueryFunctionArguments* args)
{
#ifndef ENABLE_REST_POSE_CLOSEST_POINT
    std::cout << "Rest post closest point query called with rest pose query being disabled!!!" << std::endl;
    args->query->radius = 0.f;

    return true;
#endif // DEBUG

    SurfaceMeshCloestPointQueryResult* result = (SurfaceMeshCloestPointQueryResult*)args->userPtr;
    // TetMeshFEM* pTMQuery = result->pDCD->tMeshPtrs[result->idTMQuery].get();

    VolumetricCollisionDetector* pColDet = result->pCollisionDetector;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    if (geomID != result->intersectingMeshId)
    {
        return false;
    }
    const unsigned int primID = args->primID;
    TetMeshFEM* pTMSearch = pColDet->tMeshPtrs[geomID].get();

    embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    embree::Vec3ia face(pTMSearch->surfaceFacesTetMeshVIds()(0, primID),
        pTMSearch->surfaceFacesTetMeshVIds()(1, primID), pTMSearch->surfaceFacesTetMeshVIds()(2, primID));

    embree::Vec3fa a = embree::Vec3fa::loadu(pTMSearch->mVertPos.col(face[0]).data());
    embree::Vec3fa b = embree::Vec3fa::loadu(pTMSearch->mVertPos.col(face[1]).data());
    embree::Vec3fa c = embree::Vec3fa::loadu(pTMSearch->mVertPos.col(face[2]).data());

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(queryPt, closestP);

    if (d < args->query->radius)
        // no tet traverse is needed for rest pose query
    {
        args->query->radius = d;
        result->closestFaceId = primID;
        result->faceVIds << face.x, face.y, face.z;
        result->closestFaceId = primID;
        // compute back to deformed configuration
        result->closestPt << closestP.x, closestP.y, closestP.z;

        result->barys << closestPtBarycentrics.x, closestPtBarycentrics.y, closestPtBarycentrics.z;

        result->closestPointType = pointType;

        return true; // Return true to indicate that the query radius changed.
    }

    return false;
}

GAIA::VolumetricCollisionDetector::VolumetricCollisionDetector(const VolumetricCollisionParameters& in_params)
	: params(in_params)
{
}

void tritriCollideFunc(void* userPtr, RTCCollision* collisions, unsigned int num_collisions)
{
    VolumetricCollisionDetector* pVolColDec = (VolumetricCollisionDetector*)userPtr;
    if (pVolColDec->tritriIntersectionResults.overFlow)
    {
        return;
    }
    int fId1 = collisions->primID0;
    int meshId1 = collisions->geomID0;

    int fId2 = collisions->primID1;
    int meshId2 = collisions->geomID1;
    

    TetMeshFEM::SharedPtr pTM1 = pVolColDec->tMeshPtrs[meshId1];
    TetMeshFEM::SharedPtr pTM2 = pVolColDec->tMeshPtrs[meshId2];
    
    GAIA::TriTriIntersection intersection;
    intersection.setMeshIds(meshId1, meshId2);
    const IdType* f1VIds = pTM1->surfaceFacesTetMeshVIds().col(fId1).data();
    const IdType* f2VIds = pTM2->surfaceFacesTetMeshVIds().col(fId2).data();


    bool intersectionResult = intersection.Intersect(fId1, fId2, f1VIds, f2VIds, pTM1->mVertPos, pTM2->mVertPos);

    if (intersectionResult)
    {
        TriangleCollisionResults& triCol1 = pVolColDec->triangleIntersections[meshId1][fId1];
        TriangleCollisionResults& triCol2 = pVolColDec->triangleIntersections[meshId2][fId2];
        
        int curId = pVolColDec->tritriIntersectionResults.numCollisions++;

        if (curId < pVolColDec->tritriIntersectionResults.maxNumTriTriIntersections)
        {
            pVolColDec->tritriIntersectionResults.intersections[curId] = intersection;

            int curIdTri1 = triCol1.numIntersections++;
            if (curIdTri1 < PREALLOCATED_NUM_TRI_TRI_COLLISIONS)
            {
                triCol1.intersections[curIdTri1] = curId;
            }
            else
            {
                std::cout << "[!Warning!] Triangle collision list overflow!!!!\n";
            }

            int curIdTri2 = triCol2.numIntersections++;
            if (curIdTri2 < PREALLOCATED_NUM_TRI_TRI_COLLISIONS)
            {
                triCol2.intersections[curIdTri2] = curId;
            }
            else
            {
                std::cout << "[!Warning!] Triangle collision list overflow!!!!\n";
            }
        }
        else
        {
            pVolColDec->tritriIntersectionResults.overFlow = true;
            std::cout << "[!Warning!] Tri-Tri intersection stack overflow!!!!\n";
        }
    }
    //else
    //{
    //    std::cout << "Triangles do not intersect!\n";
    //}


}


void GAIA::VolumetricCollisionDetector::initialize(std::vector<std::shared_ptr<TetMeshFEM>> tMeshes)
{
	device = rtcNewDevice(NULL);
    tMeshPtrs = tMeshes;

    surfaceMeshScene = rtcNewScene(device);
    rtcSetSceneFlags(surfaceMeshScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(surfaceMeshScene, RTC_BUILD_QUALITY_LOW);

    surfaceMeshSceneRayTracing = rtcNewScene(device);
    rtcSetSceneFlags(surfaceMeshSceneRayTracing, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
    rtcSetSceneBuildQuality(surfaceMeshSceneRayTracing, RTC_BUILD_QUALITY_LOW);

    for (int meshId = 0; meshId < tMeshes.size(); meshId++)
    {
        RTCGeometry geomRTC = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

        rtcSetSharedGeometryBuffer(geomRTC,
            RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, tMeshes[meshId]->mVertPos.data(), 0, 3 * sizeof(float),
            tMeshes[meshId]->numVertices());

        rtcSetSharedGeometryBuffer(geomRTC,
            RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, tMeshes[meshId]->surfaceFacesTetMeshVIds().data(), 0, 3 * sizeof(unsigned),
            tMeshes[meshId]->numSurfaceFaces());
        rtcSetGeometryPointQueryFunction(geomRTC, surfaceMeshClosestPointQueryFunc);

        rtcCommitGeometry(geomRTC);
        rtcAttachGeometryByID(surfaceMeshSceneRayTracing, geomRTC, meshId);
        rtcReleaseGeometry(geomRTC);


        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geom, tMeshes[meshId]->numSurfaceFaces());
        rtcSetGeometryUserData(geom, (void*)tMeshPtrs[meshId].get());
        rtcSetGeometryBoundsFunction(geom, triangle_bounds_func, nullptr);
        // rtcSetGeometryIntersectFunction(geom, triTriIntersectionFunc);

        rtcCommitGeometry(geom);
        rtcAttachGeometryByID(surfaceMeshScene, geom, meshId);
        rtcReleaseGeometry(geom);

        triangleIntersections.push_back(new TriangleCollisionResults[numSurfaceFaces(meshId)]);
        pointInclusionTestResults.emplace_back();
        pointInclusionTestResults.back().resize(tMeshPtrs[meshId]->numSurfaceVerts());

        for (size_t iSV = 0; iSV < tMeshPtrs[meshId]->numSurfaceVerts(); iSV++)
        {
            PointInclusionRTCContext& context = pointInclusionTestResults.back()[iSV];
            context.pVolCol = this;
            context.numRayIntersections.resize(tMeshPtrs.size());
        }
    }

    tritriIntersectionResults.maxNumTriTriIntersections = params.numPreallocatedTritriIntersections;
    tritriIntersectionResults.intersections.resize(params.numPreallocatedTritriIntersections);
    rtcCommitScene(surfaceMeshScene);
    rtcCommitScene(surfaceMeshSceneRayTracing);


    // initialize pointsForInclusionTest as all the points
    // for (int meshId = 0; meshId < tMeshes.size(); meshId++)
    // {
    //     pointsForInclusionTest.emplace_back();
       
    // }
}      

GAIA::VolumetricCollisionDetector::~VolumetricCollisionDetector()
{
    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        for (size_t iF = 0; iF < numSurfaceFaces(iMesh); iF++)
        {
            triangleIntersections[iMesh][iF].includedMeshes.clear();
            delete[] triangleIntersections[iMesh];
        }
    }
}

void GAIA::VolumetricCollisionDetector::updateBVH(RTCBuildQuality sceneQuality)
{
    RTCBuildQuality surfaceGeomQuality = sceneQuality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }
    rtcSetSceneBuildQuality(surfaceMeshSceneRayTracing, sceneQuality);
    rtcSetSceneBuildQuality(surfaceMeshScene, sceneQuality);


    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        RTCGeometry geomTRCSurface = rtcGetGeometry(surfaceMeshSceneRayTracing, iMesh);
        rtcSetGeometryBuildQuality(geomTRCSurface, surfaceGeomQuality);
        rtcUpdateGeometryBuffer(geomTRCSurface, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomTRCSurface);

        /*****************************************************************/

        // RTCScene surfaceMeshScene = surfaceMeshScenes[iMesh];

        // update surface Mesh
        RTCGeometry geomSurface = rtcGetGeometry(surfaceMeshScene, iMesh);
        rtcSetGeometryBuildQuality(geomSurface, surfaceGeomQuality);

        rtcUpdateGeometryBuffer(geomSurface, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomSurface);
    }
    rtcCommitScene(surfaceMeshScene);
    rtcCommitScene(surfaceMeshSceneRayTracing);

}

void GAIA::VolumetricCollisionDetector::triangleIntersectionTest()
{
    rtcCollide(surfaceMeshScene, surfaceMeshScene, tritriCollideFunc, this);
}

void GAIA::VolumetricCollisionDetector::clearTriangleIntersections()
{
    tritriIntersectionResults.clear();
    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        for (size_t iF = 0; iF < numSurfaceFaces(iMesh); iF++)
        {
            triangleIntersections[iMesh][iF].clear();
        }
    }
}

void GAIA::VolumetricCollisionDetector::findTriangleIntersectionPolygons()
{
    typedef std::pair<int, int> WorkPair;

    std::vector<WorkPair> & allWorks = tritriIntersectionResults.collidedTriangles;
    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        for (size_t iF = 0; iF < numSurfaceFaces(iMesh); iF++)
        {
            if (triangleIntersections[iMesh][iF].numIntersections)
            {
                allWorks.emplace_back(iMesh, iF);
            }
        }
    }

    cpu_parallel_for(0, allWorks.size(), [&] (int iWork) {
        clusterIntersectionPoints(allWorks[iWork].first, allWorks[iWork].second);
    });
}

void GAIA::VolumetricCollisionDetector::clusterIntersectionPoints(int meshId, int triId)
{
    TriangleCollisionResults& triColRt = triangleIntersections[meshId][triId];
    triColRt.includedVerts = 0;
    int numIntersections = std::min(triColRt.numIntersections.load(), PREALLOCATED_NUM_COLLISIONS);

    if (numIntersections == 0)
    {
        return;
    }

    
    if (numIntersections == 1)
    {
        const TriTriIntersection& tritriIntersection = tritriIntersectionResults.intersections[triColRt.intersections[0]];
        int iFace = tritriIntersection.getIFace(meshId, triId);
        for (int i = 0; i < 2; i++)
        {
            Vec3 bary;
            bary(0) = tritriIntersection.uv[iFace][i](0);
            bary(1) = tritriIntersection.uv[iFace][i](1);
            bary(2) = 1.f - bary(0) - bary(1);
            triColRt.intersectionPoints.push_back(bary);
        }
        return;
    }

    Eigen::Matrix<bool, PREALLOCATED_NUM_TRI_TRI_COLLISIONS * 2, 1> intersectionUsed = decltype(intersectionUsed)::Zero();
    Eigen::Matrix<float, 3, PREALLOCATED_NUM_TRI_TRI_COLLISIONS * 2> ptsBary;
    Eigen::Matrix<float, PREALLOCATED_NUM_TRI_TRI_COLLISIONS * 2, 3> ptsT;
    Eigen::Matrix<float, PREALLOCATED_NUM_TRI_TRI_COLLISIONS * 2, PREALLOCATED_NUM_TRI_TRI_COLLISIONS * 2> dist;

    for (int i = 0; i < numIntersections; i++)
    {
        const TriTriIntersection& tritriIntersection = tritriIntersectionResults.intersections[triColRt.intersections[i]];
        int iFace = tritriIntersection.getIFace(meshId, triId);
        ptsBary.block<2, 1>(0,i*2) = tritriIntersection.uv[iFace][0];
        ptsBary.block<2, 1>(0, i*2 + 1) = tritriIntersection.uv[iFace][1];
    }

    ptsBary.row(2) = - ptsBary.row(0) - ptsBary.row(1);
    ptsBary.row(2).array() += 1.f;

    Eigen::Matrix<float, 3, 3> triVMat;
    IdType* faceVIds = tMeshPtrs[meshId]->getSurfaceFVIdsInTetMeshVIds(triId);
    triVMat.col(0) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[0]);
    triVMat.col(1) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[1]);
    triVMat.col(2) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[2]);

    ptsT.block(0, 0, numIntersections * 2, 3) = (triVMat * ptsBary.block(0, 0, 3, numIntersections * 2)).transpose();

    for (int i = 0; i < numIntersections; i++)
    {
        dist.col(i) =
            (ptsT.rowwise() - ptsT.row(i)).matrix().rowwise().norm();
    }

    for (size_t i = 0; i < numIntersections * 2; i++)
    {
        if (!intersectionUsed(i))
        {
            triColRt.intersectionPoints.push_back(ptsBary.col(i));
        }
        for (size_t j = i+1; j < numIntersections * 2; j++)
        {
            if (dist(i, j) < params.mergeTolerance )
            {
                intersectionUsed(j) = true;
            }
        }
    }
    //// find the starting edge
    //while (true)
    //{
    //    int startEdgeId = -1;
    //    int startEdgeiFace = -1;
    //    int startEdgeEndMask = 0;
    //    int downVertMask = 7;

    //    for (int iEdge = 0; iEdge < numIntersections; ++iEdge) {

    //        const TriTriIntersection& tritriIntersection = tritriIntersectionResults.intersections[triColRt.intersections[iEdge]];
    //        int iFace = tritriIntersection.getIFace(meshId, triId);

    //        assert(iFace != -1);
    //        if ((!intersectionUsed[iEdge])
    //            && (tritriIntersection.isStartingEdge(iFace, startEdgeEndMask)))
    //        {
    //            startEdgeiFace = iFace;
    //            startEdgeId = iEdge;
    //            break;
    //        }
    //    }

    //    if (startEdgeId == -1)
    //    {
    //        std::cout << "[!Warning!] No starting edge found in a non empty intersection list!!!!\n";
    //        break;
    //    }
    //    const TriTriIntersection& tritriIntersectionStart = tritriIntersectionResults.intersections[triColRt.intersections[startEdgeId]];
    //    intersectionUsed(startEdgeId) = false;
    //    downVertMask = downVertMask & tritriIntersectionStart.outV[startEdgeiFace];
    //    // follow the edge to make a closed polygon line, which ends at another 
    //    while (true)
    //    {

    //    }
    //}
}




/* Filter callback function that gathers all hits */
void countAllHits(const struct RTCFilterFunctionNArguments* args)
{
    // std::cout << "Hit found!";

    assert(*args->valid == -1);
    PointInclusionRTCContext* context = (PointInclusionRTCContext*)args->context;

    VolumetricCollisionDetector* pVolCol = context->pVolCol;
    RTCRay* ray = (RTCRay*)args->ray;
    RTCHit* hit = (RTCHit*)args->hit;
    assert(args->N == 1);

    // The valid parameter of that structure points to an integer valid mask (0 means invalid and -1 means valid).
    args->valid[0] = 0; // ignore all hits

    /* avoid overflow of hits array */
    if (context->numHits++ >= MAX_TOTAL_HITS || (!context->sucess)) return;
    context->numHits++;

    /* check if intersection is valid */

    if (hit->u < RAY_EDGE_INTERSECTION_EPISLON || 
        hit->v < RAY_EDGE_INTERSECTION_EPISLON ||
        (1.f - hit->v - hit->u) < RAY_EDGE_INTERSECTION_EPISLON
        )
    {
        std::cout << "Hit an edge!";
        context->sucess = false;
        args->valid[0] = -1; // accept hit and redo
    }

    context->numRayIntersections(hit->geomID) += 1;

}

void GAIA::VolumetricCollisionDetector::clearpointInclusionResults()
{
    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        for (size_t iSV = 0; iSV < tMeshPtrs[iMesh]->numSurfaceVerts(); iSV++)
        {
            PointInclusionRTCContext& context = pointInclusionTestResults[iMesh][iSV];
            context.sucess = true;
            context.numHits = 0;
            context.numRayIntersections = VecDynamicI::Zero(context.numRayIntersections.rows(), 1);
        }
    }
}

void GAIA::VolumetricCollisionDetector::closestSurfacePtQuery(int iMesh, int iSV, int intersectingMeshId, SurfaceMeshCloestPointQueryResult * pResult)
{
    PointInclusionRTCContext& pointInclusionContext = pointInclusionTestResults[iMesh][iSV];
    TetMeshFEM* pTM = tMeshPtrs[iMesh].get();
    IdType vId = pTM->surfaceVIds()[iSV];
    Vec3 p = pTM->vertex(iSV);
    TetMeshFEM* pTMIntersecting = tMeshPtrs[intersectingMeshId].get();

    RTCPointQuery query;
    query.x = pTM->mVertPos(0, vId);
    query.y = pTM->mVertPos(1, vId);
    query.z = pTM->mVertPos(2, vId);
    query.radius = embree::inf;
    query.time = 0.f;

    pResult->intersectingMeshId = intersectingMeshId;
    pResult->pCollisionDetector = this;
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(surfaceMeshSceneRayTracing, &query, &context, nullptr, (void*)pResult);
    //float dis = 1e30f;

    //for (size_t iSF = 0; iSF < pTMIntersecting->numSurfaceFaces(); ++iSF) {
    //    embree::Vec3fa queryPt(p.x(), p.y(), p.z());

    //    embree::Vec3ia face(pTMIntersecting->surfaceFacesTetMeshVIds()(0, iSF),
    //        pTMIntersecting->surfaceFacesTetMeshVIds()(1, iSF), pTMIntersecting->surfaceFacesTetMeshVIds()(2, iSF));

    //    embree::Vec3fa a = embree::Vec3fa::loadu(pTMIntersecting->mVertPos.col(face[0]).data());
    //    embree::Vec3fa b = embree::Vec3fa::loadu(pTMIntersecting->mVertPos.col(face[1]).data());
    //    embree::Vec3fa c = embree::Vec3fa::loadu(pTMIntersecting->mVertPos.col(face[2]).data());

    //    ClosestPointOnTriangleType pointType;
    //    embree::Vec3fa closestPtBarycentrics;
    //    embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    //    float d = embree::distance(queryPt, closestP);

    //    if (d < dis)
    //        // no tet traverse is needed for rest pose query
    //    {
    //        dis = d;
    //        pResult->closestFaceId = iSF;
    //        pResult->faceVIds << face.x, face.y, face.z;
    //        // compute back to deformed configuration
    //        pResult->closestPt << closestP.x, closestP.y, closestP.z;

    //        pResult->barys << closestPtBarycentrics.x, closestPtBarycentrics.y, closestPtBarycentrics.z;

    //        pResult->closestPointType = pointType;

    //    }
    //}
}

void GAIA::VolumetricCollisionDetector::pointInclusionTest()
{
    
    for (size_t iMesh = 0; iMesh < tMeshPtrs.size(); iMesh++)
    {
        cpu_parallel_for(0, tMeshPtrs[iMesh]->numSurfaceVerts(), [&](int iSV) {
            PointInclusionRTCContext& context = pointInclusionTestResults[iMesh][iSV];
            rtcInitIntersectContext(&context);
            context.filter = countAllHits;

            //struct RTCRayHit
            //{
            //    struct RTCRay ray; 
            //    struct RTCHit hit;
            //};
            RTCRayHit rayhit;
            // caculate normal
            Vec3 normal;
            tMeshPtrs[iMesh]->computeVertexNormal(iSV, normal);
            rayhit.ray.dir_x = normal(0);
            rayhit.ray.dir_y = normal(1);
            rayhit.ray.dir_z = normal(2);

            int surfaceVId = tMeshPtrs[iMesh]->surfaceVIds()(iSV);
            Vec3 surfaceVert = tMeshPtrs[iMesh]->mVertPos.col(surfaceVId);
            rayhit.ray.org_x = surfaceVert(0);
            rayhit.ray.org_y = surfaceVert(1);
            rayhit.ray.org_z = surfaceVert(2);

            rayhit.ray.tnear = RAY_TRACING_NEAR;
            rayhit.ray.tfar = embree::inf;

            rtcIntersect1(surfaceMeshSceneRayTracing, &context, &rayhit);

            while (!context.sucess)
            // keep looping till success
            {
                context.sucess = true;
                context.numHits = 0;
                context.numRayIntersections = VecDynamicI::Zero(context.numRayIntersections.rows(), 1);

                Vec3 pertubatedNormal = normal + params.rayPertubation * Vec3::Random();
                pertubatedNormal /= pertubatedNormal.norm();
                rayhit.ray.dir_x = pertubatedNormal(0);
                rayhit.ray.dir_y = pertubatedNormal(1);
                rayhit.ray.dir_z = pertubatedNormal(2);

                rtcIntersect1(surfaceMeshSceneRayTracing, &context, &rayhit);
            }

            
        });
    }

}

int GAIA::VolumetricCollisionDetector::numSurfaceFaces(int meshId)
{ return tMeshPtrs[meshId]->numSurfaceFaces(); }

Vec3 GAIA::VolumetricCollisionDetector::getIntersectionPosition(int iIntersection, int triId, int meshId)
{
    Vec3 p = Vec3::Zero();

    if (iIntersection >= triangleIntersections[meshId][triId].intersectionPoints.size())
    {
        return p;
    }

    Eigen::Matrix<float, 3, 3> triVMat;
    IdType* faceVIds = tMeshPtrs[meshId]->getSurfaceFVIdsInTetMeshVIds(triId);
    triVMat.col(0) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[0]);
    triVMat.col(1) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[1]);
    triVMat.col(2) = tMeshPtrs[meshId]->mVertPos.col(faceVIds[2]);

    p = triVMat * triangleIntersections[meshId][triId].intersectionPoints[iIntersection];
     
    return p;
}

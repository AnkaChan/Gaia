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
    GAIA::TriMeshForCollision* pMesh = (GAIA::TriMeshForCollision*)args->geometryUserPtr;
    embree::BBox3fa bounds = embree::empty;

    const IdType* fVIds = pMesh->indexBuffer + 3 * args->primID;
    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->posBuffer + 3 * fVIds[0]);
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->posBuffer + 3 * fVIds[1]);
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->posBuffer + 3 * fVIds[2]);
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

    TriMeshIntersectionDetector* pColDet = result->pCollisionDetector;

    assert(args->userPtr);
    const unsigned int geomID = args->geomID;
    if (geomID != result->intersectingMeshId)
    {
        return false;
    }
    const unsigned int primID = args->primID;
    const TriMeshForCollision& mesh = pColDet->meshes[geomID];

    embree::Vec3fa queryPt(args->query->x, args->query->y, args->query->z);

    const IdType* fVIds = mesh.indexBuffer + 3 * args->primID;

    embree::Vec3fa a = embree::Vec3fa::loadu(mesh.posBuffer + 3 * fVIds[0]);
    embree::Vec3fa b = embree::Vec3fa::loadu(mesh.posBuffer + 3 * fVIds[1]);
    embree::Vec3fa c = embree::Vec3fa::loadu(mesh.posBuffer + 3 * fVIds[2]);

    ClosestPointOnTriangleType pointType;
    embree::Vec3fa closestPtBarycentrics;
    embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
    float d = embree::distance(queryPt, closestP);

    if (d < args->query->radius)
        // no tet traverse is needed for rest pose query
    {
        args->query->radius = d;
        result->closestFaceId = primID;
        result->faceVIds << fVIds[0], fVIds[1], fVIds[2];
        result->closestFaceId = primID;
        // compute back to deformed configuration
        result->closestPt << closestP.x, closestP.y, closestP.z;

        result->barys << closestPtBarycentrics.x, closestPtBarycentrics.y, closestPtBarycentrics.z;

        result->closestPointType = pointType;

        return true; // Return true to indicate that the query radius changed.
    }

    return false;
}

GAIA::TriMeshIntersectionDetector::TriMeshIntersectionDetector(const TriMeshIntersectionDetectorParameters& in_params)
	: params(in_params)
{
}

void tritriCollideFunc(void* userPtr, RTCCollision* collisions, unsigned int num_collisions)
{
    TriMeshIntersectionDetector* pVolColDec = (TriMeshIntersectionDetector*)userPtr;
    if (pVolColDec->tritriIntersectionResults.overFlow)
    {
        return;
    }
    int fId1 = collisions->primID0;
    int meshId1 = collisions->geomID0;

    int fId2 = collisions->primID1;
    int meshId2 = collisions->geomID1;
    

    const TriMeshForCollision & mesh1 = pVolColDec->meshes[meshId1];
    const TriMeshForCollision & mesh2 = pVolColDec->meshes[meshId2];
    
    GAIA::TriTriIntersection intersection;
    intersection.setMeshIds(meshId1, meshId2);
    const IdType* f1VIds = mesh1.getFaceVIds(fId1);
    const IdType* f2VIds = mesh2.getFaceVIds(fId2);


    bool intersectionResult = intersection.Intersect(fId1, fId2, f1VIds, f2VIds, mesh1.posBuffer, mesh2.posBuffer);

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


void GAIA::TriMeshIntersectionDetector::initialize(std::vector<TriMeshForCollision> meshes)
{
	device = rtcNewDevice(NULL);
    meshes = std::move(meshes);

    triMeshIntersectionScene = rtcNewScene(device);
    rtcSetSceneFlags(triMeshIntersectionScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(triMeshIntersectionScene, RTC_BUILD_QUALITY_LOW);

    size_t numAllFace = 0;

    for (int meshId = 0; meshId < meshes.size(); meshId++)
    {
        numAllFace += meshes[meshId].numFaces;

        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geom, meshes[meshId].numFaces);
        rtcSetGeometryUserData(geom, (void*)&meshes[meshId]);
        rtcSetGeometryBoundsFunction(geom, triangle_bounds_func, nullptr);
        // rtcSetGeometryIntersectFunction(geom, triTriIntersectionFunc);

        rtcCommitGeometry(geom);
        rtcAttachGeometryByID(triMeshIntersectionScene, geom, meshId);
        rtcReleaseGeometry(geom);

        triangleIntersections.push_back(new TriangleCollisionResults[numSurfaceFaces(meshId)]);
    }

    tritriIntersectionResults.maxNumTriTriIntersections = params.tritriIntersectionsPreallocatedRatio* numAllFace;
    tritriIntersectionResults.intersections.resize(tritriIntersectionResults.maxNumTriTriIntersections);
    rtcCommitScene(triMeshIntersectionScene);

    if (params.enbalePointInclusionTest)
    {
        triMeshSceneRayTracing = rtcNewScene(device);
        rtcSetSceneFlags(triMeshSceneRayTracing, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST | RTC_SCENE_FLAG_CONTEXT_FILTER_FUNCTION);
        rtcSetSceneBuildQuality(triMeshSceneRayTracing, RTC_BUILD_QUALITY_LOW);
        for (int meshId = 0; meshId < meshes.size(); meshId++)
        {
            RTCGeometry geomRTC = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_TRIANGLE);

            rtcSetSharedGeometryBuffer(geomRTC,
                RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, meshes[meshId].posBuffer, 0, 3 * sizeof(float),
                meshes[meshId].numVertices);

            rtcSetSharedGeometryBuffer(geomRTC,
                RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, meshes[meshId].indexBuffer, 0, 3 * sizeof(unsigned),
                meshes[meshId].numFaces);
            rtcSetGeometryPointQueryFunction(geomRTC, surfaceMeshClosestPointQueryFunc);

            rtcCommitGeometry(geomRTC);
            rtcAttachGeometryByID(triMeshSceneRayTracing, geomRTC, meshId);
            rtcReleaseGeometry(geomRTC);

            pointInclusionTestResults.emplace_back();
            pointInclusionTestResults.back().resize(meshes[meshId].numFaces);

            for (size_t iSV = 0; iSV < meshes[meshId].numVertices; iSV++)
            {
                PointInclusionRTCContext& context = pointInclusionTestResults.back()[iSV];
                context.pVolCol = this;
                context.numRayIntersections.resize(meshes.size());
            }
        }

        rtcCommitScene(triMeshSceneRayTracing);
    }


    // initialize pointsForInclusionTest as all the points
    // for (int meshId = 0; meshId < tMeshes.size(); meshId++)
    // {
    //     pointsForInclusionTest.emplace_back();
       
    // }
}      

GAIA::TriMeshIntersectionDetector::~TriMeshIntersectionDetector()
{
    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
    {
        for (size_t iF = 0; iF < numSurfaceFaces(iMesh); iF++)
        {
            triangleIntersections[iMesh][iF].includedMeshes.clear();
            delete[] triangleIntersections[iMesh];
        }
    }
}

void GAIA::TriMeshIntersectionDetector::updateBVH(RTCBuildQuality sceneQuality)
{
    RTCBuildQuality surfaceGeomQuality = sceneQuality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }
    rtcSetSceneBuildQuality(triMeshIntersectionScene, sceneQuality);

    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
    {

        // update surface Mesh
        RTCGeometry geomSurface = rtcGetGeometry(triMeshIntersectionScene, iMesh);
        rtcSetGeometryBuildQuality(geomSurface, surfaceGeomQuality);

        rtcUpdateGeometryBuffer(geomSurface, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomSurface);
    }
    rtcCommitScene(triMeshIntersectionScene);

    if (params.enbalePointInclusionTest)
    {
        rtcSetSceneBuildQuality(triMeshSceneRayTracing, sceneQuality);
        for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
        {
            RTCGeometry geomTRCSurface = rtcGetGeometry(triMeshSceneRayTracing, iMesh);
            rtcSetGeometryBuildQuality(geomTRCSurface, surfaceGeomQuality);
            rtcUpdateGeometryBuffer(geomTRCSurface, RTC_BUFFER_TYPE_VERTEX, 0);
            rtcCommitGeometry(geomTRCSurface);
        }
        rtcCommitScene(triMeshSceneRayTracing);
    }
}

void GAIA::TriMeshIntersectionDetector::triangleIntersectionTest()
{
    rtcCollide(triMeshIntersectionScene, triMeshIntersectionScene, tritriCollideFunc, this);
}

void GAIA::TriMeshIntersectionDetector::clearTriangleIntersections()
{
    tritriIntersectionResults.clear();
    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
    {
        for (size_t iF = 0; iF < numSurfaceFaces(iMesh); iF++)
        {
            triangleIntersections[iMesh][iF].clear();
        }
    }
}

void GAIA::TriMeshIntersectionDetector::findTriangleIntersectionPolygons()
{
    typedef std::pair<int, int> WorkPair;

    std::vector<WorkPair> & allWorks = tritriIntersectionResults.collidedTriangles;
    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
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

void GAIA::TriMeshIntersectionDetector::clusterIntersectionPoints(int meshId, int triId)
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
    IdType* faceVIds = meshes[meshId].getFaceVIds(triId);
    triVMat.col(0) << meshes[meshId].getVertex(faceVIds[0])[0], meshes[meshId].getVertex(faceVIds[0])[1], meshes[meshId].getVertex(faceVIds[0])[2];
    triVMat.col(1) << meshes[meshId].getVertex(faceVIds[1])[0], meshes[meshId].getVertex(faceVIds[1])[1], meshes[meshId].getVertex(faceVIds[1])[2];
    triVMat.col(2) << meshes[meshId].getVertex(faceVIds[2])[0], meshes[meshId].getVertex(faceVIds[2])[1], meshes[meshId].getVertex(faceVIds[2])[2];

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

    TriMeshIntersectionDetector* pVolCol = context->pVolCol;
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

void GAIA::TriMeshIntersectionDetector::clearpointInclusionResults()
{
    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
    {
        for (size_t iSV = 0; iSV < meshes[iMesh].numVertices; iSV++)
        {
            PointInclusionRTCContext& context = pointInclusionTestResults[iMesh][iSV];
            context.sucess = true;
            context.numHits = 0;
            context.numRayIntersections = VecDynamicI::Zero(context.numRayIntersections.rows(), 1);
        }
    }
}

void GAIA::TriMeshIntersectionDetector::closestSurfacePtQuery(int iMesh, int vId, int intersectingMeshId, SurfaceMeshCloestPointQueryResult * pResult)
{
    PointInclusionRTCContext& pointInclusionContext = pointInclusionTestResults[iMesh][vId];
    const TriMeshForCollision& mesh = meshes[iMesh];
    Vec3 p = mesh.getVertexVec3(vId);

    RTCPointQuery query;
    query.x = p(0);
    query.y = p(1);
    query.z = p(2);
    query.radius = embree::inf;
    query.time = 0.f;

    pResult->intersectingMeshId = intersectingMeshId;
    pResult->pCollisionDetector = this;
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);
    rtcPointQuery(triMeshSceneRayTracing, &query, &context, nullptr, (void*)pResult);
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

void GAIA::TriMeshIntersectionDetector::pointInclusionTest()
{
    
    for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
    {
        cpu_parallel_for(0, meshes[iMesh].numVertices, [&](int iSV) {
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

            std::cout << "Normal computation is not implemented!!!";
            RELEASE_ASSERT(false);

            //meshes[iMesh]->computeVertexNormal(iSV, normal);
            //rayhit.ray.dir_x = normal(0);
            //rayhit.ray.dir_y = normal(1);
            //rayhit.ray.dir_z = normal(2);

            //int surfaceVId = tMeshPtrs[iMesh]->surfaceVIds()(iSV);
            //Vec3 surfaceVert = tMeshPtrs[iMesh]->mVertPos.col(surfaceVId);
            //rayhit.ray.org_x = surfaceVert(0);
            //rayhit.ray.org_y = surfaceVert(1);
            //rayhit.ray.org_z = surfaceVert(2);

            //rayhit.ray.tnear = RAY_TRACING_NEAR;
            //rayhit.ray.tfar = embree::inf;

            //rtcIntersect1(triMeshSceneRayTracing, &context, &rayhit);

            //while (!context.sucess)
            //// keep looping till success
            //{
            //    context.sucess = true;
            //    context.numHits = 0;
            //    context.numRayIntersections = VecDynamicI::Zero(context.numRayIntersections.rows(), 1);

            //    Vec3 pertubatedNormal = normal + params.rayPertubation * Vec3::Random();
            //    pertubatedNormal /= pertubatedNormal.norm();
            //    rayhit.ray.dir_x = pertubatedNormal(0);
            //    rayhit.ray.dir_y = pertubatedNormal(1);
            //    rayhit.ray.dir_z = pertubatedNormal(2);

            //    rtcIntersect1(triMeshSceneRayTracing, &context, &rayhit);
            //}
            
        });
    }

}

int GAIA::TriMeshIntersectionDetector::numSurfaceFaces(int meshId)
{ return meshes[meshId].numFaces; }

Vec3 GAIA::TriMeshIntersectionDetector::getIntersectionPosition(int iIntersection, int triId, int meshId)
{
    Vec3 p = Vec3::Zero();

    if (iIntersection >= triangleIntersections[meshId][triId].intersectionPoints.size())
    {
        return p;
    }

    Eigen::Matrix<float, 3, 3> triVMat;
    IdType* faceVIds = meshes[meshId].getFaceVIds(triId);
    triVMat.col(0) = meshes[meshId].getVertexVec3(faceVIds[0]);
    triVMat.col(1) = meshes[meshId].getVertexVec3(faceVIds[1]);
    triVMat.col(2) = meshes[meshId].getVertexVec3(faceVIds[2]);

    p = triVMat * triangleIntersections[meshId][triId].intersectionPoints[iIntersection];
     
    return p;
}

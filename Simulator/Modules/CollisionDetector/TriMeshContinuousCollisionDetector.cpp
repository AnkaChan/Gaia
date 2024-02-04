#include "TriMeshContinuousCollisionDetector.h"
#include "../TriMesh/TriMesh.h"

#include "CCDSolver.h"

#define DEBUG_CCD
typedef double CCDDType;

using namespace GAIA;


template<typename DType>
inline void setVert(cy::Vec3<DType>& v, const Vec3& v_in) {
    v.x = (CCDDType)v_in.x();
    v.y = (CCDDType)v_in.y();
    v.z = (CCDDType)v_in.z();
}


GAIA::TriMeshContinuousCollisionDetector::TriMeshContinuousCollisionDetector(const CollisionDetectionParamters& in_params) 
	: params(in_params)
{
}


void movingFaceBoundsFuncTriMesh(const struct RTCBoundsFunctionArguments* args)
{
    const CCDGeometry* pGeom = (CCDGeometry*)args->geometryUserPtr;

    const GAIA::TriMeshFEM* pMesh = pGeom->pMesh;
    embree::BBox3fa bounds = embree::empty;
    const int faceId = args->primID;
    const GAIA::IdType* face = pMesh->facePos.col(faceId).data();


    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->vertex(face[0]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->vertex(face[1]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->vertex(face[2]).data());
    bounds.extend(a);
    
    // use extenal prev pos intead of the pMesh->positionsPrev
    const TVerticesMat* pPrevPos = pGeom->pPrevPos;
    a = embree::Vec3fa::loadu(pPrevPos->col(face[0]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pPrevPos->col(face[1]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pPrevPos->col(face[2]).data());
    bounds.extend(a);

    *(embree::BBox3fa*)args->bounds_o = bounds;
}

void movingVertsBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const CCDGeometry* pGeom = (CCDGeometry*)args->geometryUserPtr;

    const GAIA::TriMeshFEM* pMesh = pGeom->pMesh;
    embree::BBox3fa bounds = embree::empty;
    const int vId = args->primID;
    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->vertex(vId).data());
    bounds.extend(a);

    // use extenal prev pos intead of the pMesh->positionsPrev
    const TVerticesMat* pPrevPos = pGeom->pPrevPos;
    a = embree::Vec3fa::loadu(pPrevPos->col(vId).data());
    bounds.extend(a);

    *(embree::BBox3fa*)args->bounds_o = bounds;
}

void movingEdgesBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const CCDGeometry* pGeom = (CCDGeometry*)args->geometryUserPtr;

    const GAIA::TriMeshFEM* pMesh = pGeom->pMesh;
    embree::BBox3fa bounds = embree::empty;
    const int eId = args->primID;

    const EdgeInfo & eInfo = pMesh->pTopology->edgeInfos[eId];
    embree::Vec3fa a = embree::Vec3fa::loadu(pMesh->vertex(eInfo.eV1).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pMesh->vertex(eInfo.eV2).data());
    bounds.extend(a);

    // use extenal prev pos intead of the pMesh->positionsPrev
    const TVerticesMat* pPrevPos = pGeom->pPrevPos;
    a = embree::Vec3fa::loadu(pPrevPos->col(eInfo.eV1).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pPrevPos->col(eInfo.eV2).data());
    bounds.extend(a);

    //if (eId == 129 || eId == 217 || eId == 303 || eId == 216)
    //if (eId == 1440 || eId == 1351)
    //{
    //    std::cout << "Edge " << eId << " bounds: " << bounds.lower << " " << bounds.upper << std::endl;
    //}

    *(embree::BBox3fa*)args->bounds_o = bounds;
}



//bool continuousTriPointIntersectionFunc(RTCPointQueryFunctionArguments* args)
//{
//    GAIA::CollisionDetectionResult* result = (GAIA::CollisionDetectionResult*)args->userPtr;
//    unsigned int  primID = args->primID;
//    unsigned int  geomID = args->geomID;
//
//    int intersectedMeshId = geomID;
//    GAIA::TriMeshContinuousCollisionDetector* pCCD = (GAIA::TriMeshContinuousCollisionDetector*)result->pDetector;
//
//    GAIA::TriMeshFEM* pMIntersected = pCCD->meshPtrs[geomID].get();
//    GAIA::TriMeshFEM* pMQuery = pCCD->meshPtrs[result->idTMQuery].get();
//    GAIA::IdType* face = pMIntersected->facePos.col(primID).data();
//
//    if (intersectedMeshId == result->idTMQuery)
//    {
//        // to do detect if the vertex is on the fIntersected
//        // if it does, this face will not be counted as intersected 
//        for (int iFV = 0; iFV <3; iFV++)
//        {
//            if (face[iFV] == result->idVQuery)
//            {
//                return false;
//            }
//        }
//    }
//
//    cy::Vec3<CCDDType> fvs[2][3];
//
//    // face from intersected mesh
//    int iV = 0;
//    for (int iFV = 0; iFV < 3; iFV++)
//    {
//        fvs[0][iV].x = (CCDDType)pMIntersected->positionsPrev(0, face[iFV]);
//        fvs[0][iV].y = (CCDDType)pMIntersected->positionsPrev(1, face[iFV]);
//        fvs[0][iV].z = (CCDDType)pMIntersected->positionsPrev(2, face[iFV]);
//
//        fvs[1][iV].x = (CCDDType)pMIntersected->positions()(0, face[iFV]);
//        fvs[1][iV].y = (CCDDType)pMIntersected->positions()(1, face[iFV]);
//        fvs[1][iV].z = (CCDDType)pMIntersected->positions()(2, face[iFV]);
//        ++iV;
//    }
//
//    // point from query mesh
//    cy::Vec3<CCDDType> p[2];
//
//    p[0].x = pMQuery->positionsPrev(0, result->idVQuery);
//    p[0].y = pMQuery->positionsPrev(1, result->idVQuery);
//    p[0].z = pMQuery->positionsPrev(2, result->idVQuery);
//
//    p[1].x = pMQuery->positions()(0, result->idVQuery);
//    p[1].y = pMQuery->positions()(1, result->idVQuery);
//    p[1].z = pMQuery->positions()(2, result->idVQuery);
//
//    // IntersectContinuousTriPoint(float& tt, Vec3f const x[2][3], Vec3f const p[2])
//
//    double tt = -1;
//    cy::Vec3<CCDDType> barycentrics;
//    if (cy::IntersectContinuousTriPoint(tt, fvs, p, barycentrics)) {
//        float penetrationDepth = tt * (p[1] - p[0]).Length();
//        if (!result->numIntersections()) {
//            result->closestSurfaceFaceId.push_back(primID);
//            // result->intersectedTets.push_back(-1);
//            result->intersectedTMeshIds.push_back(intersectedMeshId);
//            result->shortestPathFound.push_back(true);
//            result->penetrationDepth = penetrationDepth;
//        }
//        else if(result->penetrationDepth > penetrationDepth)
//        {
//            result->closestSurfaceFaceId[0] = (primID);
//            result->closestSurfaceFaceId[0] = (primID);
//            // result->intersectedTets.push_back(-1);
//            result->intersectedTMeshIds[0] = (intersectedMeshId);
//            result->shortestPathFound[0] = (true);
//            result->penetrationDepth = penetrationDepth;
//        }
//        // result->closestSurfacePtBarycentrics.push_back({ (float)barycentrics.x, (float)barycentrics.y, (float)barycentrics.z });
//        // result->intersectedFaces.push_back(&fIntersected);
//
//#ifdef DEBUG_CCD
//        //std::cout << "Collision detected for Vectex: " << result->pVQuery->pTetMeshV->id()
//        //    << " from tetmesh " << result->pMQuery->pTetMesh->id
//        //    << " at time step: " << tt << std::endl;
//        CPoint penetratePointFace(0, 0, 0);
//
//        iV = 0;
//        for (M::VPtr pV : It::FVIterator(&fIntersected))
//        {
//            CPoint& preVPos = pV->prevPos;
//            CPoint& endVPos = pV->point();
//            CPoint fvPosIntersection = endVPos * tt + preVPos * (1 - tt);
//
//            penetratePointFace = penetratePointFace + fvPosIntersection * barycentricsP[iV];
//            ++iV;
//        }
//
//
//        CPoint  penetratePoint = result->pVQuery->prevPos * (1 - tt)
//            + result->pVQuery->point() * tt;
//
//        //std::cout << "Face point when intersection: " << penetratePointFace;
//        //std::cout << "| Vertex position when intersection: " << penetratePoint << std::endl;
//        static double worst = 0;
//        double dis = (penetratePoint - penetratePointFace).norm();
//        if (dis > 1e-6)
//        {
//            if (worst < dis)
//            {
//                worst = dis;
//            }
//            std::cout << "Distance: " << (penetratePoint - penetratePointFace).norm() << " | worst: " << worst << "\n";
//            std::cout << "In comparison, length of point trajectory: " << (result->pVQuery->prevPos - result->pVQuery->point()).norm()
//                << " | face edge length: \n";
//
//            for (M::EPtr pE : It::FEIterator(&fIntersected))
//            {
//                std::cout << M::edgeLength(pE) << " ";
//
//            }
//            std::cout << "\n";
//            assert(false);
//
//        }
//
//
//#endif // DEBUG_CCD
//
//        // the return valud means whether the search radius have changed
//        return false;
//    }
//    else
//    {
//        return false;
//    }
//}

void GAIA::TriMeshContinuousCollisionDetector::initialize(std::vector<std::shared_ptr<TriMeshFEM>> meshes, std::vector<TVerticesMat*> prevPoses)
{
    // add all the tet mesh to a single scene for collision detection
    device = rtcNewDevice(NULL);

    // create a scene for the moving faces
    triangleTrajectoriesScene = rtcNewScene(device);
    rtcSetSceneFlags(triangleTrajectoriesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(triangleTrajectoriesScene, RTC_BUILD_QUALITY_LOW);

    // create a scene for the moving vertices
    vertexTrajectoriesScene = rtcNewScene(device);
    rtcSetSceneFlags(vertexTrajectoriesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(vertexTrajectoriesScene, RTC_BUILD_QUALITY_LOW);

    // create a scene for the moving edges
    edgeTrajectoriesScene = rtcNewScene(device);
    rtcSetSceneFlags(edgeTrajectoriesScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(edgeTrajectoriesScene, RTC_BUILD_QUALITY_LOW);

    //mMeshPrevPosHandles.resize(meshPtrs.size());

    assert(prevPoses.size() == meshes.size());

    size_t numVerts = 0, numEdges = 0;

    for (int meshId = 0; meshId < meshes.size(); meshId++)
    {
        ccdGeometries.emplace_back();
        TriMeshFEM* pMesh = meshes[meshId].get();
        ccdGeometries.back().pMesh = pMesh;
        ccdGeometries.back().pPrevPos = prevPoses[meshId];

        CCDGeometry* pGeom = &ccdGeometries.back();
        unsigned int geomId = meshId;

        /* triangles */
        RTCGeometry geomTris = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geomTris, pMesh->numFaces());
        rtcSetGeometryBoundsFunction(geomTris, movingFaceBoundsFuncTriMesh, nullptr);

        rtcSetGeometryUserData(geomTris, (void*)(pGeom));
        // no need to set up query function, because we are colliding the geometry with a moving point scene
        rtcCommitGeometry(geomTris);
        rtcAttachGeometryByID(triangleTrajectoriesScene, geomTris, geomId);
        rtcReleaseGeometry(geomTris);

        /* vertices */
        RTCGeometry geomVerts = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);

        rtcSetGeometryUserPrimitiveCount(geomVerts, pMesh->numVertices());
        rtcSetGeometryBoundsFunction(geomVerts, movingVertsBoundsFunc, nullptr);
        rtcSetGeometryUserData(geomVerts, (void*)(pGeom));
        rtcCommitGeometry(geomVerts);
        rtcAttachGeometryByID(vertexTrajectoriesScene, geomVerts, geomId);
        rtcReleaseGeometry(geomVerts);
        numVerts += pMesh->numVertices();

        /* edges */
        RTCGeometry geomEdges = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geomEdges, pMesh->numEdges());
        rtcSetGeometryBoundsFunction(geomEdges, movingEdgesBoundsFunc, nullptr);
        rtcSetGeometryUserData(geomEdges, (void*)(pGeom));
        rtcCommitGeometry(geomEdges);
        rtcAttachGeometryByID(edgeTrajectoriesScene, geomEdges, geomId);
        rtcReleaseGeometry(geomEdges);
        numEdges += pMesh->numEdges();
    }
    rtcCommitScene(triangleTrajectoriesScene);
    rtcCommitScene(vertexTrajectoriesScene);
    rtcCommitScene(edgeTrajectoriesScene);

    vfCollisionResults.initialize(preAllocationRatio * numVerts);
    eeCollisionResults.initialize(preAllocationRatio * numEdges);
    
}


void GAIA::TriMeshContinuousCollisionDetector::updateBVH(RTCBuildQuality quality)
{
    RTCBuildQuality sceneQuality = quality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }

    rtcSetSceneBuildQuality(triangleTrajectoriesScene, sceneQuality);
    rtcSetSceneBuildQuality(vertexTrajectoriesScene, sceneQuality);
    rtcSetSceneBuildQuality(edgeTrajectoriesScene, sceneQuality);

    for (int meshId = 0; meshId < ccdGeometries.size(); meshId++)
    {
        // get the tet geom buffer
        // unsign        // unsigned int geoId = mMeshGeometryIds[meshId];
        int geoId = meshId;

        RTCGeometry geomTri = rtcGetGeometry(triangleTrajectoriesScene, geoId);
        rtcSetGeometryBuildQuality(geomTri, quality);
        rtcUpdateGeometryBuffer(geomTri, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomTri);

        RTCGeometry geomVert = rtcGetGeometry(vertexTrajectoriesScene, geoId);
        rtcSetGeometryBuildQuality(geomVert, quality);
        rtcUpdateGeometryBuffer(geomVert, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomVert);

        RTCGeometry geomEdge = rtcGetGeometry(edgeTrajectoriesScene, geoId);
        rtcSetGeometryBuildQuality(geomEdge, quality);
        rtcUpdateGeometryBuffer(geomEdge, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geomEdge);
    }

    rtcCommitScene(triangleTrajectoriesScene);
    rtcCommitScene(vertexTrajectoriesScene);
    rtcCommitScene(edgeTrajectoriesScene);
}

void vfCollideFunc(void* userPtr, RTCCollision* collisions, unsigned int num_collisions)
{
    TriMeshContinuousCollisionDetector* pCCD = (TriMeshContinuousCollisionDetector*)userPtr;
    if (pCCD->vfCollisionResults.overflow)
    {
        return;
    }
    const int vId1 = collisions->primID0;
    const int meshId1 = collisions->geomID0;

    const int fId2 = collisions->primID1;
    const int meshId2 = collisions->geomID1;

    const TriMeshFEM::Ptr pMesh1 = pCCD->ccdGeometries[meshId1].pMesh;
    const TVerticesMat* pPrevPos1 = pCCD->ccdGeometries[meshId1].pPrevPos;
    const TriMeshFEM::Ptr pMesh2 = pCCD->ccdGeometries[meshId2].pMesh;
    const TVerticesMat* pPrevPos2 = pCCD->ccdGeometries[meshId2].pPrevPos;
    const GAIA::IdType* face = pMesh2->facePos.col(fId2).data();

    // filter out the self collision
    if (meshId1 == meshId2)
    {
        for (int iFV = 0; iFV < 3; iFV++)
        {
            if (vId1 == face[iFV])
            {
                return;
            }
        }
    }

    // face from mesh 2
    cy::Vec3<CCDDType> fvs[2][3];
    for (int iFV = 0; iFV < 3; iFV++)
    {
        const Vec3 prevPos = pPrevPos2->col(face[iFV]);
        fvs[0][iFV].x = (CCDDType)(prevPos.x());
        fvs[0][iFV].y = (CCDDType)(prevPos.y());
        fvs[0][iFV].z = (CCDDType)(prevPos.z());

        const Vec3 pos = pMesh2->vertex(face[iFV]);
        fvs[1][iFV].x = (CCDDType)(pos.x());
        fvs[1][iFV].y = (CCDDType)(pos.y());
        fvs[1][iFV].z = (CCDDType)(pos.z());
    }

    // point from mesh 1
    cy::Vec3<CCDDType> p[2];
    const Vec3 vPrevPos = pPrevPos1->col(vId1);
    const Vec3 vPos = pMesh1->vertex(vId1);

    p[0].x = (CCDDType)vPrevPos.x();
    p[0].y = (CCDDType)vPrevPos.y();
    p[0].z = (CCDDType)vPrevPos.z();

    p[1].x = (CCDDType)vPos.x();
    p[1].y = (CCDDType)vPos.y();
    p[1].z = (CCDDType)vPos.z();


    cy::Vec3<CCDDType> barycentrics;
    CCDDType tt = -1;
    if (cy::IntersectContinuousTriPoint(tt, fvs, p, barycentrics))
    {
        const int curId = pCCD->vfCollisionResults.numCollisions++;
        if (curId < pCCD->vfCollisionResults.maxNumTriTriIntersections)
        {
            VFCollision& curCollision = pCCD->vfCollisionResults.collisions[curId];
            curCollision.t = tt;
            curCollision.barycentrics(0) = (FloatingType)barycentrics.x;
            curCollision.barycentrics(1) = (FloatingType)barycentrics.y;
            curCollision.barycentrics(2) = (FloatingType)barycentrics.z;

            curCollision.faceMeshId = meshId2;
            curCollision.faceId = fId2;

            curCollision.vertexMeshId = meshId1;
            curCollision.vertexId = vId1;

#ifdef RECORD_COLLIDING_POINT
            //Vec3 c(0, 0, 0);
            //for (int iFV = 0; iFV < 3; iFV++)
            //{
            //    const Vec3 prevPos = pPrevPos2->col(face[iFV]);
            //    c += (1. - tt) *  prevPos * barycentrics[iFV];

            //    const Vec3 pos = pMesh2->vertex(face[iFV]);
            //    c += tt * pos * barycentrics[iFV];
            //}

            Vec3 c = (1. - tt) * vPrevPos + tt * vPos;
            curCollision.c = c;
#endif //RECORD_COLLIDING_POINT
        }
		else
		{
			pCCD->vfCollisionResults.overflow = true;
			return;
        }


#ifdef DEBUG_CCD
        int iV = 0;
        Vec3 penetratePointFacePrev(0, 0, 0);
        Vec3 penetratePointFaceNow(0, 0, 0);
        for (int iFV = 0; iFV < 3; iFV++)
        {
            const Vec3 prevPosF = pPrevPos2->col(face[iFV]);
            penetratePointFacePrev += prevPosF * barycentrics[iV];

            const Vec3 posF = pMesh2->vertex(face[iFV]);
            penetratePointFaceNow += posF * barycentrics[iV];

            ++iV;
        }

        Vec3 penetratePointFace = (1. - tt) * penetratePointFacePrev + tt * penetratePointFaceNow;

        Vec3 penetratePoint = (1.-tt) * vPrevPos + tt * vPos;

        static double worst = 0;
        double dis = (penetratePoint - penetratePointFace).norm();
        if (dis > 1e-3f)
        {
            if (worst < dis)
            {
                worst = dis;
            }
            std::cout << "Error! Inaccurate VF CCD encountered! Distance: " << (penetratePoint - penetratePointFace).norm() << " | worst: " << worst << "\n";
            //std::cout << "In comparison, length of point trajectory: " << (result->pVQuery->prevPos - result->pVQuery->point()).norm()
            //    << " | face edge length: \n";

            //for (M::EPtr pE : It::FEIterator(&fIntersected))
            //{
            //    std::cout << M::edgeLength(pE) << " ";

            //}
            //std::cout << "\n";
            assert(false);
        }
#endif // DEBUG_CCD
    }
}

void eeCollideFunc(void* userPtr, RTCCollision* collisions, unsigned int num_collisions)
{
    TriMeshContinuousCollisionDetector* pCCD = (TriMeshContinuousCollisionDetector*)userPtr;
    if (pCCD->eeCollisionResults.overflow)
    {
        return;
    }

    const int eId1 = collisions->primID0;
    const int meshId1 = collisions->geomID0;

    const int eId2 = collisions->primID1;
    const int meshId2 = collisions->geomID1;

    const TriMeshFEM::Ptr pMesh1 = pCCD->ccdGeometries[meshId1].pMesh;
    const TVerticesMat* pPrevPos1 = pCCD->ccdGeometries[meshId1].pPrevPos;
    const TriMeshFEM::Ptr pMesh2 = pCCD->ccdGeometries[meshId2].pMesh;
    const TVerticesMat* pPrevPos2 = pCCD->ccdGeometries[meshId2].pPrevPos;

    const EdgeInfo & edgeInfo1 = pMesh1->pTopology->edgeInfos[eId1];
    const EdgeInfo & edgeInfo2 = pMesh2->pTopology->edgeInfos[eId2];

    //if (eId1 == 1440 && eId2 == 1351)
    //{
    //    std::cout << "Edge " << eId1 << " and " << eId2 << " collide!\n";
    //}

    //if (eId1 == 1351 && eId2 == 1440)
    //{
    //    std::cout << "Edge " << eId1 << " and " << eId2 << " collide!\n";
    //}

    // filter out the self collision
    if (meshId1 == meshId2)
    {
        // same edge
        if (eId1 == eId2)
        {
            return;
        }
        // adjacent edges
        if (edgeInfo1.eV1 == edgeInfo2.eV1
            || edgeInfo1.eV1 == edgeInfo2.eV2
            || edgeInfo1.eV2 == edgeInfo2.eV1
            || edgeInfo1.eV2 == edgeInfo2.eV2)
        {
            return;
        }
    }


    
    cy::Vec3<CCDDType> e0[2][2];
    cy::Vec3<CCDDType> e1[2][2];

    // edge from mesh 1
    const Vec3 e0V0Prev = pPrevPos1->col(edgeInfo1.eV1);
    setVert(e0[0][0], e0V0Prev);

    const Vec3 e0V1Prev = pPrevPos1->col(edgeInfo1.eV2);
    setVert(e0[0][1], e0V1Prev);

    const Vec3 e0V0 = pMesh1->vertex(edgeInfo1.eV1);
    setVert(e0[1][0], e0V0);

    const Vec3 e0V1 = pMesh1->vertex(edgeInfo1.eV2);
    setVert(e0[1][1], e0V1);

    // edge from mesh 2
    const Vec3 e1V0Prev = pPrevPos2->col(edgeInfo2.eV1);
    setVert(e1[0][0], e1V0Prev);

    const Vec3 e1V1Prev = pPrevPos2->col(edgeInfo2.eV2);
    setVert(e1[0][1], e1V1Prev);

    const Vec3 e1V0 = pMesh2->vertex(edgeInfo2.eV1);
    setVert(e1[1][0], e1V0);

    const Vec3 e1V1 = pMesh2->vertex(edgeInfo2.eV2);
	setVert(e1[1][1], e1V1);


    CCDDType tt = -1;

    if (cy::IntersectContinuousEdgeEdge(tt, e0, e1))
    {
        const int curId = pCCD->eeCollisionResults.numCollisions++;
        if (curId < pCCD->eeCollisionResults.maxNumTriTriIntersections)
        {
            const Vec3 e0V0_collide = (1. - tt) * e0V0Prev + tt * e0V0;
            const Vec3 e0V1_collide = (1. - tt) * e0V1Prev + tt * e0V1;

            const Vec3 e1V0_collide = (1. - tt) * e1V0Prev + tt * e1V0;
            const Vec3 e1V1_collide = (1. - tt) * e1V1Prev + tt * e1V1;

            const embree::Vec3fa p1 = embree::Vec3fa::loadu(e0V0_collide.data());
            const embree::Vec3fa p2 = embree::Vec3fa::loadu(e0V1_collide.data());

            const embree::Vec3fa q1 = embree::Vec3fa::loadu(e1V0_collide.data());
            const embree::Vec3fa q2 = embree::Vec3fa::loadu(e1V1_collide.data());
            embree::Vec3fa c1, c2;
            FloatingType mua, mub;
            get_closest_points_between_segments(p1, p2, q1, q2, c1, c2, mua, mub);

            EECollision& curCollision = pCCD->eeCollisionResults.collisions[curId];
            curCollision.miu1 = mua;
            curCollision.miu2 = mub;

            curCollision.t = tt;

#ifdef RECORD_COLLIDING_POINT
            Vec3 c;
            c << c1.x, c1.y, c1.z;
            curCollision.c = c;
#endif //RECORD_COLLIDING_POINT

#ifdef DEBUG_CCD
            static double worst = 0;
            double dis = embree::length(c1 - c2);
            if (dis > 1e-4)
            {
                if (worst < dis)
                {
                    worst = dis;
                }
                std::cout << "Error! Inaccurate EE CCD encountered! Distance: " 
                    << dis << " | worst: " << worst 
                    << " | t: " << tt << " | mua: " << mua << " | mub: " << mub 
                    << " between: " << eId1 << " and " << eId2 << "\n";



                //std::cout << "In comparison, length of point trajectory: " << (result->pVQuery->prevPos - result->pVQuery->point()).norm()
                //    << " | face edge length: \n";

                //for (M::EPtr pE : It::FEIterator(&fIntersected))
                //{
                //    std::cout << M::edgeLength(pE) << " ";

                //}
                //std::cout << "\n";
                //assert(false);
            }
        }
        else
        {
            pCCD->eeCollisionResults.overflow = true;
            return;
        }


        
#endif // DEBUG_CCD
    }
}

bool GAIA::TriMeshContinuousCollisionDetector::continuousCollisionDetection()
{
    do
    {
        vfCollisionResults.clear();
        rtcCollide(vertexTrajectoriesScene, triangleTrajectoriesScene, vfCollideFunc, this);
    } while (vfCollisionResults.overflow);

    do
    {
        eeCollisionResults.clear();
        rtcCollide(edgeTrajectoriesScene, edgeTrajectoriesScene, eeCollideFunc, this);
    } while (eeCollisionResults.overflow);

    //rtcCollide
    return false;
}


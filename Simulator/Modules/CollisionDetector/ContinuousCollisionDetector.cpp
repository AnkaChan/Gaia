#include "ContinuousCollisionDetector.h"
#include "../TetMesh/TetMeshFEM.h"

#include "CCDSolver.h"

#include "CollisionGeometry.h"

//typedef double CCDDType;
typedef GAIA::FloatingType CCDDType;

// #define DEBUG_CCD

GAIA::ContinuousCollisionDetector::ContinuousCollisionDetector(const CollisionDetectionParamters& in_params) 
	: params(in_params)
{
}


void movingFaceBoundsFunc(const struct RTCBoundsFunctionArguments* args)
{
    const GAIA::TetMeshFEM* pTM = (const GAIA::TetMeshFEM*)args->geometryUserPtr;

    //Eigen::Matrix<GAIA::FloatingType, 6, 3> faceTrajPrism;

    //for (int iV = 0; iV < 3; iV++)
    //{
    //    faceTrajPrism.row(iV) = pTM->mVertPos.col(pTM->surfaceFacesTetMeshVIds()(iV, args->primID)).transpose();
    //}

    //for (int iV = 0; iV < 3; iV++)
    //{
    //    faceTrajPrism.row(iV+3) = pTM->mVertPrevPos.col(pTM->surfaceFacesTetMeshVIds()(iV, args->primID)).transpose();
    //}

    //RTCBounds* bounds_o = args->bounds_o;

    //bounds_o->upper_x = faceTrajPrism.col(0).maxCoeff();
    //bounds_o->upper_y = faceTrajPrism.col(1).maxCoeff();
    //bounds_o->upper_z = faceTrajPrism.col(2).maxCoeff();
    //bounds_o->lower_x = faceTrajPrism.col(0).minCoeff();
    //bounds_o->lower_y = faceTrajPrism.col(1).minCoeff();
    //bounds_o->lower_z = faceTrajPrism.col(2).minCoeff();


    embree::BBox3fa bounds = embree::empty;
    const int faceId = args->primID;
    const GAIA::IdType* face = pTM->surfaceFacesTetMeshVIds().col(args->primID).data();


    embree::Vec3fa a = embree::Vec3fa::loadu(pTM->vertex(face[0]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pTM->vertex(face[1]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(pTM->vertex(face[2]).data());
    bounds.extend(a);

    // use extenal prev pos intead of the pMesh->positionsPrev
    const GAIA::TVerticesMat & prevPos = pTM->mVertPrevPos;
    a = embree::Vec3fa::loadu(prevPos.col(face[0]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(prevPos.col(face[1]).data());
    bounds.extend(a);
    a = embree::Vec3fa::loadu(prevPos.col(face[2]).data());
    bounds.extend(a);

    *(embree::BBox3fa*)args->bounds_o = bounds;
}


bool continuousTriPointIntersectionFunc(RTCPointQueryFunctionArguments* args)
{
    GAIA::CollisionDetectionResult* result = (GAIA::CollisionDetectionResult*)args->userPtr;
    unsigned int  primID = args->primID;
    unsigned int  geomID = args->geomID;

    int intersectedMeshId = geomID;
    GAIA::ContinuousCollisionDetector* pCCD = (GAIA::ContinuousCollisionDetector*)result->pDetector;

    GAIA::TetMeshFEM* pMIntersected = pCCD->tMeshPtrs[geomID].get();
    GAIA::TetMeshFEM* pMQuery = pCCD->tMeshPtrs[result->idTMQuery].get();
    const GAIA::IdType* face = pMIntersected->surfaceFacesTetMeshVIds().col(primID).data();

    if (intersectedMeshId == result->idTMQuery)
    {
        // to do detect if the vertex is on the fIntersected
        // if it does, this face will not be counted as intersected 
        for (int iFV = 0; iFV <3; iFV++)
        {
            if (face[iFV] == result->idVQuery)
            {
                return false;
            }
        }
    }

    cy::Vec3<CCDDType> fvs[2][3];

    // face from intersected mesh
    int iV = 0;
    for (int iFV = 0; iFV < 3; iFV++)
    {
        fvs[0][iV].x = (CCDDType)pMIntersected->mVertPrevPos(0, face[iFV]);
        fvs[0][iV].y = (CCDDType)pMIntersected->mVertPrevPos(1, face[iFV]);
        fvs[0][iV].z = (CCDDType)pMIntersected->mVertPrevPos(2, face[iFV]);

        fvs[1][iV].x = (CCDDType)pMIntersected->mVertPos(0, face[iFV]);
        fvs[1][iV].y = (CCDDType)pMIntersected->mVertPos(1, face[iFV]);
        fvs[1][iV].z = (CCDDType)pMIntersected->mVertPos(2, face[iFV]);
        ++iV;
    }

    // point from query mesh
    cy::Vec3<CCDDType> p[2];

    p[0].x = pMQuery->mVertPrevPos(0, result->idVQuery);
    p[0].y = pMQuery->mVertPrevPos(1, result->idVQuery);
    p[0].z = pMQuery->mVertPrevPos(2, result->idVQuery);

    p[1].x = pMQuery->mVertPos(0, result->idVQuery);
    p[1].y = pMQuery->mVertPos(1, result->idVQuery);
    p[1].z = pMQuery->mVertPos(2, result->idVQuery);

    // IntersectContinuousTriPoint(float& tt, Vec3f const x[2][3], Vec3f const p[2])

    CCDDType tt = -1;
    cy::Vec3<CCDDType> barycentrics;
    if (cy::IntersectContinuousTriPoint<CCDDType>(tt, fvs, p, barycentrics)) {
        GAIA::FloatingType penetrationDepth = tt * (p[1] - p[0]).Length();

        bool needUpdate = false;
        if (!result->numIntersections()) {
            result->collidingPts.emplace_back();
            needUpdate = true;
        }
        else if (result->penetrationDepth > penetrationDepth)
        {
            needUpdate = true;
        }
        GAIA::CollidingPointInfo& colldingPt = result->collidingPts.back();

        if (needUpdate)
        {
                GAIA::Vec3 penetratePoint = (1. - tt) * pMQuery->vertexPrevPos(result->idVQuery)
                + tt * pMQuery->vertex(result->idVQuery);

                colldingPt.closestSurfaceFaceId = primID;
                // result->intersectedTets.push_back(-1);
                colldingPt.intersectedMeshId = intersectedMeshId;
                colldingPt.shortestPathFound = true;
                colldingPt.closestSurfacePtBarycentrics <<
                    (GAIA::FloatingType)barycentrics[0], (GAIA::FloatingType)barycentrics[1], (GAIA::FloatingType)barycentrics[2];

                colldingPt.closestSurfacePt = penetratePoint;
                colldingPt.closestPointType = GAIA::ClosestPointOnTriangleType::AtInterior;

                result->penetrationDepth = penetrationDepth;

                GAIA::Vec3 contactNormal;
                if (pCCD->params.computeContactNormal)
                {
                    GAIA::computeContactNormalTetMesh(*result, 0, contactNormal, pCCD->tMeshPtrs);
                }
                else
                {
                    contactNormal << 0.f, 0.f, 0.f;
                }
                colldingPt.closestPointNormal = contactNormal;
        }

        //    CollidingPointInfo& colldingPt = result->collidingPts.back();

        //    result->closestSurfaceFaceId.push_back(primID);
        //    // result->intersectedTets.push_back(-1);
        //    result->intersectedTMeshIds.push_back(intersectedMeshId);
        //    result->shortestPathFound.push_back(true);
        //    result->closestSurfacePtBarycentrics.push_back(
        //        { (GAIA::FloatingType)barycentrics[0] , (GAIA::FloatingType)barycentrics[1], (GAIA::FloatingType)barycentrics[2] }
        //    );
        //    result->closestSurfacePts.push_back(penetratePoint);
        //    result->closestPointType.push_back(GAIA::ClosestPointOnTriangleType::AtInterior);

        //    result->penetrationDepth = penetrationDepth;

        //    GAIA::Vec3 contactNormal;
        //    if (pCCD->params.computeContactNormal)
        //    {
        //        GAIA::computeContactNormalTetMesh(*result, 0, contactNormal, pCCD->tMeshPtrs);
        //    }
        //    else
        //    {
        //        contactNormal << 0.f, 0.f, 0.f;
        //    }
        //    result->closestPointNormals.push_back(contactNormal);
        //}
        //else if(result->penetrationDepth > penetrationDepth)
        //{
        //    result->closestSurfaceFaceId[0] = (primID);
        //    // result->intersectedTets.push_back(-1);
        //    result->intersectedTMeshIds[0] = (intersectedMeshId);
        //    result->closestSurfacePtBarycentrics[0] 
        //        << (GAIA::FloatingType)barycentrics[0], (GAIA::FloatingType)barycentrics[1], (GAIA::FloatingType)barycentrics[2];
        //    result->closestSurfacePts[0] = penetratePoint;
        //    result->shortestPathFound[0] = true;
        //    result->closestPointType[0] = GAIA::ClosestPointOnTriangleType::AtInterior;
        //    result->penetrationDepth = penetrationDepth;

        //    GAIA::Vec3 contactNormal;
        //    if (pCCD->params.computeContactNormal)
        //    {
        //        GAIA::computeContactNormalTetMesh(*result, 0, contactNormal, pCCD->tMeshPtrs);
        //    }
        //    else
        //    {
        //        contactNormal << 0.f, 0.f, 0.f;
        //    }
        //    result->closestPointNormals[0] = contactNormal;
        //}

        // result->closestSurfacePtBarycentrics.push_back({ (float)barycentrics.x, (float)barycentrics.y, (float)barycentrics.z });
        // result->intersectedFaces.push_back(&fIntersected);

#ifdef DEBUG_CCD
        GAIA::Vec3 penetratePoint = (1. - tt) * pMQuery->vertexPrevPos(result->idVQuery)

        GAIA::Vec3 penetratePointFacePrev(0, 0, 0);
        GAIA::Vec3 penetratePointFaceNow(0, 0, 0);
        int iV = 0;
        for (int iFV = 0; iFV < 3; iFV++)
        {
            const GAIA::Vec3 prevPosF = pMIntersected->mVertPrevPos.col(face[iFV]);
            penetratePointFacePrev += prevPosF * barycentrics[iV];

            const GAIA::Vec3 posF = pMIntersected->vertex(face[iFV]);
            penetratePointFaceNow += posF * barycentrics[iV];

            ++iV;
        }
        GAIA::Vec3 penetratePointFace = (1. - tt) * penetratePointFacePrev + tt * penetratePointFaceNow;

        static double worst = 0;
        double dis = (penetratePoint - penetratePointFace).norm();
        // std::cout << "CCD collision distance: " << dis << "\n";
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

        // the return valud means whether the search radius have changed
        return false;
    }
    else
    {
        return false;
    }


}

void GAIA::ContinuousCollisionDetector::initialize(std::vector<std::shared_ptr<TetMeshFEM>> tMeshes)
{
    numFaces = 0;
    tMeshPtrs = tMeshes;
    // add all the tet mesh to a single scene for collision detection
    device = rtcNewDevice(NULL);
    surfaceTriangleTrajectoryScene = rtcNewScene(device);
    rtcSetSceneFlags(surfaceTriangleTrajectoryScene, RTC_SCENE_FLAG_DYNAMIC | RTC_SCENE_FLAG_ROBUST);
    rtcSetSceneBuildQuality(surfaceTriangleTrajectoryScene, RTC_BUILD_QUALITY_LOW);

    //mMeshPrevPosHandles.resize(meshPtrs.size());

    for (int meshId = 0; meshId < tMeshes.size(); meshId++)
    {
        // add the dynmanic prop handle for previous frame position
        //pM->addVProp(mMeshPrevPosHandles[tMId]);

        TetMeshFEM* pTM = tMeshes[meshId].get();
        /* Uses custom geometry */
        RTCGeometry geom = rtcNewGeometry(device, RTC_GEOMETRY_TYPE_USER);
        rtcSetGeometryUserPrimitiveCount(geom, pTM->numSurfaceFaces());
        rtcSetGeometryBoundsFunction(geom, movingFaceBoundsFunc, nullptr);

        rtcSetGeometryUserData(geom, (void*)(pTM));

        rtcSetGeometryPointQueryFunction(geom, continuousTriPointIntersectionFunc);
        rtcCommitGeometry(geom);
        unsigned int geomId = meshId;
        rtcAttachGeometryByID(surfaceTriangleTrajectoryScene, geom, geomId);
        //mMeshGeoIdToPointerMap[geomId] = pM;
        //mMeshPtrToMeshIndexMap[pM] = tMId;
        rtcReleaseGeometry(geom);

    }
    rtcCommitScene(surfaceTriangleTrajectoryScene);

    
}


void GAIA::ContinuousCollisionDetector::updateBVH(RTCBuildQuality quality)
{
    RTCBuildQuality sceneQuality = quality;
    if (sceneQuality == RTC_BUILD_QUALITY_REFIT) {
        sceneQuality = RTC_BUILD_QUALITY_LOW;
    }

    rtcSetSceneBuildQuality(surfaceTriangleTrajectoryScene, sceneQuality);

    for (int meshId = 0; meshId < tMeshPtrs.size(); meshId++)
    {
        // get the tet geom buffer
        // unsign        // unsigned int geoId = mMeshGeometryIds[meshId];
        int geoId = meshId;
        TetMeshFEM* pTM = tMeshPtrs[meshId].get();

        RTCGeometry geom = rtcGetGeometry(surfaceTriangleTrajectoryScene, geoId);

        if (tMeshPtrs[meshId]->activeForCollision) {
            rtcEnableGeometry(geom);
        }
        else
        {
            rtcDisableGeometry(geom);
            continue;
        }
        rtcSetGeometryBuildQuality(geom, quality);
        rtcUpdateGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0);
        rtcCommitGeometry(geom);

    }

    rtcCommitScene(surfaceTriangleTrajectoryScene);
}


bool GAIA::ContinuousCollisionDetector::vertexContinuousCollisionDetection(int32_t vId, int32_t tMeshId, CollisionDetectionResult* pResult)
{
    RTCPointQuery query;
    RTCPointQueryContext context;
    rtcInitPointQueryContext(&context);

    pResult->clear();
    pResult->idVQuery = vId;
    pResult->idTMQuery = tMeshId;
    pResult->pDetector = (void*)this;
    pResult->fromCCD = true;

    //query.radius = (prevPosVQuery - pVQuery->point()).norm();
    // radius should be half of the length that the vertex moved
    TetMeshFEM* pTM = tMeshPtrs[tMeshId].get();
    Vec3 middleTrajectory = 0.5f * (pTM->mVertPos.col(vId) + pTM->mVertPrevPos.col(vId));

    query.radius = (pTM->mVertPos.col(vId) - pTM->mVertPrevPos.col(vId)).norm() * 0.5f;
    query.time = 0.f;

    query.x = middleTrajectory(0);
    query.y = middleTrajectory(1);
    query.z = middleTrajectory(2);


    rtcPointQuery(surfaceTriangleTrajectoryScene, &query, &context, nullptr, (void*)pResult);

    //rtcCollide
    return false;
}


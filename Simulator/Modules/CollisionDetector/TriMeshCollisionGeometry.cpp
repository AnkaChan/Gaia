#pragma once
#include "TriMeshCollisionGeometry.h"

#define ABSOLUTE_RELAXIATION (1e-6f)


bool GAIA::checkFeasibleRegion(const embree::Vec3fa& p, const TriMeshFEM* pTM, int32_t faceId, ClosestPointOnTriangleType pointType, float feasibleRegionEpsilon)
{
    bool inFeasibleRegion = true;
    const int32_t* faceVIds = pTM->facePos.col(faceId).data();
    switch (pointType)
    {
    case ClosestPointOnTriangleType::AtInterior:
        // this is automatically satisfied
        inFeasibleRegion = true;
        break;
    case ClosestPointOnTriangleType::AtAB:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 0, feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtBC:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 1, feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtAC:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 2, feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtA:
        inFeasibleRegion = checkVertexFeasibleRegion(p, pTM, faceVIds[0], feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtB:
        inFeasibleRegion = checkVertexFeasibleRegion(p, pTM, faceVIds[1], feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtC:
        inFeasibleRegion = checkVertexFeasibleRegion(p, pTM, faceVIds[2], feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::NotFound:
        inFeasibleRegion = false;
        break;
    default:
        break;
    }
    return inFeasibleRegion;
}

//bool GAIA::checkEdgeFeasibleRegion(const embree::Vec3fa& p, const TriMeshFEM* pMesh, int32_t faceId, int32_t edgeId,
//    int32_t edgeVId1, int32_t edgeVId2, float feasibleRegionEpsilon)
//{
//    int32_t neighborFaceId = pMesh->pTopology->faces3NeighborFaces(edgeId, faceId);
//
//    embree::Vec3fa v1 = loadVertexPos(pMesh, edgeVId1);
//    embree::Vec3fa v2 = loadVertexPos(pMesh, edgeVId2);
// 
//    embree::Vec3fa AP = p - v1;
//    embree::Vec3fa BP = p - v2;
//    embree::Vec3fa AB = v2 - v1;
//    embree::Vec3fa BA = -AB;
//    
//    float relaxed = -feasibleRegionEpsilon - ABSOLUTE_RELAXIATION;
//
//    if (embree::dot(AP, AB) < relaxed) {
//        return false;
//    }
//    //CPoint BP = p - B;
//
//    if (embree::dot(BP, BA) < relaxed) {
//        return false;
//    }
//    
//    //embree::Vec3fa fNormal1 = faceNormal(pMesh, faceId);
//
//    //embree::Vec3fa nAB = cross(AB, fNormal1);
//    //if (embree::dot(AP, nAB) < relaxed) {
//    //    return false;
//    //}
//
//
//
//    if (neighborFaceId != -1) {
//		// non-boundary edge, need to check the other side
//        
//        embree::Vec3fa fNormal2 = faceNormal(pMesh, neighborFaceId);
//
//        embree::Vec3fa nBA = cross(BA, fNormal2);
//        if (embree::dot(AP, nBA) < relaxed) {
//            return false;
//        }
//    }
//
//    return true;
//}

bool GAIA::checkEdgeFeasibleRegion(const embree::Vec3fa& p, const TriMeshFEM* pMesh, int32_t faceId,
    int32_t edgeOrderInFace, float feasibleRegionEpsilon)
{
    IdType edgeId = pMesh->pTopology->faces3NeighborEdges(edgeOrderInFace, faceId);
    const EdgeInfo& edgeInfo = pMesh->pTopology->edgeInfos[edgeId];
    CIdType edgeVId1 = edgeInfo.eV1;
    CIdType edgeVId2 = edgeInfo.eV2;

    embree::Vec3fa v1 = loadVertexPos(pMesh, edgeVId1);
    embree::Vec3fa v2 = loadVertexPos(pMesh, edgeVId2);

    embree::Vec3fa AP = p - v1;
    embree::Vec3fa BP = p - v2;
    embree::Vec3fa AB = v2 - v1;
    float relaxed = feasibleRegionEpsilon + ABSOLUTE_RELAXIATION;

    if (embree::dot(AP, AB) < -relaxed) {
        return false;
    }

    if (embree::dot(BP, AB) > relaxed) {
        return false;
    }

    CFloatingType ABSqrNorm = embree::sqr_length(AB);
    if (ABSqrNorm < CMP_EPSILON2)
    {
        return false;
    }

    if (edgeInfo.eV12Next != -1)
    {
        embree::Vec3fa v12Next = loadVertexPos(pMesh, edgeInfo.eV12Next);
        CFloatingType t = embree::dot(AB, v12Next - v1) / ABSqrNorm;

        embree::Vec3fa perpendiculaFoot = v1 + t * AB; 
        // RELEASE_ASSERT(abs(embree::dot(v12Next - perpendiculaFoot, AB)) < CMP_EPSILON);

        // the normal of the divider plane is (v12Next - perpendiculaFoot).normalized()
        // to be located in the feasible region p has to be located in the different side of divider plane as v12Next
        if (embree::dot(v12Next - perpendiculaFoot, p - perpendiculaFoot) > relaxed)
        {
            return false;
        }
    }

    if (edgeInfo.eV21Next != -1)
    {
        embree::Vec3fa v21Next = loadVertexPos(pMesh, edgeInfo.eV21Next);
        CFloatingType t = embree::dot(AB, v21Next - v1) / ABSqrNorm;

        embree::Vec3fa perpendiculaFoot = v1 + t * AB;
        // RELEASE_ASSERT(abs(embree::dot(v21Next - perpendiculaFoot, AB)) < CMP_EPSILON);
        // the normal of the divider plane is (v21Next - perpendiculaFoot).normalized()
        // to be located in the feasible region p has to be located in the different side of divider plane as v21Next
        if (embree::dot(v21Next - perpendiculaFoot, p - perpendiculaFoot) > relaxed)
        {
            return false;
        }
    }

    return true;
}

bool GAIA::checkVertexFeasibleRegion(const embree::Vec3fa& p, const TriMeshFEM* pMesh, int32_t vId, float feasibleRegionEpsilon)
{
    embree::Vec3fa A = loadVertexPos(pMesh, vId);;

    embree::Vec3fa AP = p - A;

    size_t numNeiVertices = pMesh->numNeiVertices(vId);

    for (int32_t iVNei = 0; iVNei < numNeiVertices; ++iVNei) {
        size_t neiVId = pMesh->getVertexIthNeiVertex(vId, iVNei);
        embree::Vec3fa B = loadVertexPos(pMesh, neiVId);;

        embree::Vec3fa BA = A - B;

        //// skip feasible region filtering for inverted surface parts
        //if (!pTetM->DCDEnabled(pVNei->pTetMeshV))
        //{
        //    return true;
        //}

        float relaxed = -feasibleRegionEpsilon - ABSOLUTE_RELAXIATION;

        // AP  * BA > 0 means P is on the demanded side of plane passing through A whose normal is BA 
        // add a little margin to make the determination more conservative: 
        // instead of AP * BA <= 0 we use AP * BA <= epsilon, where epsilon < 0
        if (embree::dot(AP, BA) < relaxed) {
            return false;
        }
    }

    return true;
}



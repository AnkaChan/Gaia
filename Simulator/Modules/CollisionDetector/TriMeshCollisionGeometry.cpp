#pragma once
#include "TriMeshCollisionGeometry.h"

#define ABSOLUTE_RELAXIATION (1e-6f)


bool GAIA::checkFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pTM, int32_t faceId, ClosestPointOnTriangleType pointType, float feasibleRegionEpsilon)
{
    bool inFeasibleRegion = true;
    int32_t* faceVIds = pTM->facePos.col(faceId).data();
    switch (pointType)
    {
    case ClosestPointOnTriangleType::AtInterior:
        // this is automatically satisfied
        inFeasibleRegion = true;
        break;
    case ClosestPointOnTriangleType::AtAB:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 0, faceVIds[0], faceVIds[1], feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtBC:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 1, faceVIds[1], faceVIds[2], feasibleRegionEpsilon);
        break;
    case ClosestPointOnTriangleType::AtAC:
        inFeasibleRegion = checkEdgeFeasibleRegion(p, pTM, faceId, 2, faceVIds[2], faceVIds[0], feasibleRegionEpsilon);
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

bool GAIA::checkEdgeFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pMesh, int32_t faceId, int32_t edgeId,
    int32_t edgeVId1, int32_t edgeVId2, float feasibleRegionEpsilon)
{
    int32_t neighborFaceId = pMesh->pTopology->faces3NeighborFaces(edgeId, faceId);

    embree::Vec3fa v1 = loadVertexPos(pMesh, edgeVId1);
    embree::Vec3fa v2 = loadVertexPos(pMesh, edgeVId2);
 
    embree::Vec3fa AP = p - v1;
    embree::Vec3fa BP = p - v2;
    embree::Vec3fa AB = v2 - v1;
    embree::Vec3fa BA = -AB;
    
    float relaxed = -feasibleRegionEpsilon - ABSOLUTE_RELAXIATION;

    if (embree::dot(AP, AB) < relaxed) {
        return false;
    }
    //CPoint BP = p - B;

    if (embree::dot(BP, BA) < relaxed) {
        return false;
    }
    
    embree::Vec3fa fNormal1 = faceNormal(pMesh, faceId);

    embree::Vec3fa nAB = cross(AB, fNormal1);
    if (embree::dot(AP, nAB) < relaxed) {
        return false;
    }

    if (neighborFaceId != -1) {
		// non-boundary edge, need to check the other side
        
        embree::Vec3fa fNormal2 = faceNormal(pMesh, neighborFaceId);

        embree::Vec3fa nBA = cross(BA, fNormal2);
        if (embree::dot(AP, nBA) < relaxed) {
            return false;
        }
    }

    return true;
}

bool GAIA::checkVertexFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pMesh, int32_t vId, float feasibleRegionEpsilon)
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



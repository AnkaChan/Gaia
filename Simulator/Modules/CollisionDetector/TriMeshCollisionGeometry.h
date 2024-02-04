#pragma once
#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "CollisionDetertionParameters.h"

#include "../TriMesh/TriMesh.h"

namespace GAIA {

    inline embree::Vec3fa loadVertexPos(TriMeshFEM* pMesh, int32_t vId)
    {
        return embree::Vec3fa::loadu(pMesh->vertex(vId).data());
    }

    embree::Vec3fa faceNormal(TriMeshFEM* pMesh, int32_t faceId);
    embree::Vec3fa faceOrientedArea(TriMeshFEM* pMesh, int32_t faceId);


    bool checkFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pMesh, int32_t faceId,
        ClosestPointOnTriangleType pointType, float feasibleRegionEpsilon);

    bool checkEdgeFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pMesh, int32_t faceId,
        int32_t edgeId, int32_t edgeVId1, int32_t edgeVId2, float feasibleRegionEpsilon);

    bool checkVertexFeasibleRegion(const embree::Vec3fa& p, TriMeshFEM* pTM, int32_t vId, float feasibleRegionEpsilon);


    inline embree::Vec3fa faceOrientedArea(TriMeshFEM* pMesh, int32_t faceId)
    {
        embree::Vec3fa a = loadVertexPos(pMesh, pMesh->facePosVId(faceId, 0)),
            b = loadVertexPos(pMesh, pMesh->facePosVId(faceId, 1)),
            c = loadVertexPos(pMesh, pMesh->facePosVId(faceId, 2));

        embree::Vec3fa ab = b - a;
        embree::Vec3fa ac = c - a;

        embree::Vec3fa orientedArea = embree::cross(ab, ac);

        return orientedArea;
    }

    inline embree::Vec3fa GAIA::faceNormal(TriMeshFEM* pMesh, int32_t faceId)
    {
        embree::Vec3fa normal = faceOrientedArea(pMesh, faceId);
        normal = embree::normalize(normal);

        return normal;
    }
}
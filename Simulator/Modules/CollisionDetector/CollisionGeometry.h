#pragma once
#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "CollisionDetertionParameters.h"

#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"

#define d_of(m, n, o, p) ((m.x - n.x) * (o.x - p.x) + (m.y - n.y) * (o.y - p.y) + (m.z - n.z) * (o.z - p.z))


namespace GAIA {

	inline bool pointInTet(embree::Vec3fa const& p, embree::Vec3fa const& a, embree::Vec3fa const& b, embree::Vec3fa const& c, embree::Vec3fa const& d) {
		const embree::Vec3fa AB = b - a;
		const embree::Vec3fa AC = c - a;
		const embree::Vec3fa AD = d - a;

		FloatingType tetOrientedVol = dot(AB, cross(AC, AD));
		
		if (abs(tetOrientedVol) < CMP_EPSILON)
		{
			return false;
		}

		const embree::Vec3fa* vs[4] = {
			&a, &b, &c, &d
		};

		const int32_t order[4][3] = { { 1, 2, 3 },{ 2, 0, 3 },{ 0, 1, 3 },{ 1, 0, 2 } };

		for (int32_t i = 0; i < 4; ++i) {
			embree::Vec3fa const& v1 = *(vs[order[i][0]]);
			embree::Vec3fa const& v2 = *(vs[order[i][1]]);
			embree::Vec3fa const& v3 = *(vs[order[i][2]]);

			const embree::Vec3fa v = p - v1;

			if (dot(v, cross(v2 - v1, v3 - v2)) * tetOrientedVol >= 0)
			{
				return false;
			}
		}

		return true;
	}

	inline embree::Vec3fa closestPointTriangle(embree::Vec3fa const& p, embree::Vec3fa const& a, embree::Vec3fa const& b,
		embree::Vec3fa const& c, embree::Vec3fa& baryCentrics, ClosestPointOnTriangleType& pointType)
    {
        const embree::Vec3fa ab = b - a;
        const embree::Vec3fa ac = c - a;
        const embree::Vec3fa ap = p - a;

        const float d1 = dot(ab, ap);
        const float d2 = dot(ac, ap);
        if (d1 <= 0.f && d2 <= 0.f) {
            pointType = ClosestPointOnTriangleType::AtA;
            baryCentrics = embree::Vec3fa(1.f, 0.f, 0.f);

            return a;
        }

        const embree::Vec3fa bp = p - b;
        const float d3 = dot(ab, bp);
        const float d4 = dot(ac, bp);
        if (d3 >= 0.f && d4 <= d3) {
            pointType = ClosestPointOnTriangleType::AtB;
            baryCentrics = embree::Vec3fa(0.f, 1.f, 0.f);
            return b;
        }

        const embree::Vec3fa cp = p - c;
        const float d5 = dot(ab, cp);
        const float d6 = dot(ac, cp);
        if (d6 >= 0.f && d5 <= d6) {
            pointType = ClosestPointOnTriangleType::AtC;
            baryCentrics = embree::Vec3fa(0.f, 0.f, 1.f);
            return c;
        }

        const float vc = d1 * d4 - d3 * d2;
        if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f)
        {
            const float v = d1 / (d1 - d3);
            pointType = ClosestPointOnTriangleType::AtAB;
            baryCentrics = embree::Vec3fa(1.0f - v, v, 0.f);
            return a + v * ab;
        }

        const float vb = d5 * d2 - d1 * d6;
        if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f)
        {
            const float v = d2 / (d2 - d6);
            pointType = ClosestPointOnTriangleType::AtAC;
            baryCentrics = embree::Vec3fa(1.0f - v, 0.f, v);
            return a + v * ac;
        }

        const float va = d3 * d6 - d5 * d4;
        if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f)
        {
            pointType = ClosestPointOnTriangleType::AtBC;
            const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));
            baryCentrics = embree::Vec3fa(0.f, 1.f - v, v);
            return b + v * (c - b);
        }

        const float denom = 1.f / (va + vb + vc);
        const float v = vb * denom;
        const float w = vc * denom;
        pointType = ClosestPointOnTriangleType::AtInterior;
        baryCentrics = embree::Vec3fa(1.f - v - w, v, w);
        return a + v * ab + w * ac;
    }


	inline void get_closest_points_between_segments(const embree::Vec3fa& p1, const embree::Vec3fa& p2, const embree::Vec3fa& q1, 
		const embree::Vec3fa& q2, embree::Vec3fa& c1, embree::Vec3fa& c2, FloatingType& mua, FloatingType& mub) {
		// Calculate the parametric position on the 2 curves, mua and mub.
		FloatingType denominator1 = (d_of(p2, p1, p2, p1) * d_of(q2, q1, q2, q1) - d_of(q2, q1, p2, p1) * d_of(q2, q1, p2, p1));
		FloatingType denominator2 = d_of(q2, q1, q2, q1);
		if (Utility::isZeroApprox(denominator1) || Utility::isZeroApprox(denominator2)) {
			// The lines are parallel.
			mua = 0.f;
			mub = 0.f;
			c1 = p1;
			c2 = q1;
		}
		else
		{
			mua = (d_of(p1, q1, q2, q1) * d_of(q2, q1, p2, p1) - d_of(p1, q1, p2, p1) * d_of(q2, q1, q2, q1))
				/ denominator1;
			mub = (d_of(p1, q1, q2, q1) + mua * d_of(q2, q1, p2, p1)) / denominator2;

			// Clip the value between [0..1] constraining the solution to lie on the original curves.
			if (mua < 0) {
				mua = 0;
			}
			if (mub < 0) {
				mub = 0;
			}
			if (mua > 1) {
				mua = 1;
			}
			if (mub > 1) {
				mub = 1;
			}
			c1 = lerp(p1, p2, mua);
			c2 = lerp(q1, q2, mub);
		}

	}

	inline FloatingType get_closest_distance_between_segments(const embree::Vec3fa& p_from_a, const embree::Vec3fa& p_to_a, const embree::Vec3fa& p_from_b, const embree::Vec3fa& p_to_b) {
		embree::Vec3fa u = p_to_a - p_from_a;
		embree::Vec3fa v = p_to_b - p_from_b;
		embree::Vec3fa w = p_from_a - p_to_a;
		FloatingType a = embree::dot(u, u);; // u.dot(u); // Always >= 0
		FloatingType b = embree::dot(u, v);; // u.dot(v);
		FloatingType c = embree::dot(v, v);; // v.dot(v); // Always >= 0
		FloatingType d = embree::dot(u, w);; // u.dot(w);
		FloatingType e = embree::dot(v, w);; // v.dot(w);
		FloatingType D = a * c - b * b; // Always >= 0
		FloatingType sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
		FloatingType tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0

		// Compute the line parameters of the two closest points.
		if (D < CMP_EPSILON) { // The lines are almost parallel.
			sN = 0.0; // Force using point P0 on segment S1
			sD = 1.0; // to prevent possible division by 0.0 later.
			tN = e;
			tD = c;
		}
		else { // Get the closest points on the infinite lines
			sN = (b * e - c * d);
			tN = (a * e - b * d);
			if (sN < 0.0) { // sc < 0 => the s=0 edge is visible.
				sN = 0.0;
				tN = e;
				tD = c;
			}
			else if (sN > sD) { // sc > 1 => the s=1 edge is visible.
				sN = sD;
				tN = e + b;
				tD = c;
			}
		}

		if (tN < 0.0) { // tc < 0 => the t=0 edge is visible.
			tN = 0.0;
			// Recompute sc for this edge.
			if (-d < 0.0) {
				sN = 0.0;
			}
			else if (-d > a) {
				sN = sD;
			}
			else {
				sN = -d;
				sD = a;
			}
		}
		else if (tN > tD) { // tc > 1 => the t=1 edge is visible.
			tN = tD;
			// Recompute sc for this edge.
			if ((-d + b) < 0.0) {
				sN = 0;
			}
			else if ((-d + b) > a) {
				sN = sD;
			}
			else {
				sN = (-d + b);
				sD = a;
			}
		}

		// Finally do the division to get sc and tc.
		sc = (Utility::isZeroApprox(sN) ? 0.0 : sN / sD);
		tc = (Utility::isZeroApprox(tN) ? 0.0 : tN / tD);

		// Get the difference of the two closest points.
		embree::Vec3fa dP = w + (sc * u) - (tc * v); // = S1(sc) - S2(tc)

		return embree::length(dP); // Return the closest distance.
	}

	inline void computeContactNormalTetMesh(CollisionDetectionResult& colResult, int32_t iIntersection, 
		Vec3& normal, std::vector<std::shared_ptr<TetMeshFEM>>& tMeshPtrs)
	{
		CollidingPointInfo& collidingPt = colResult.collidingPts[iIntersection];
		ClosestPointOnTriangleType pointType = collidingPt.closestPointType;
		int32_t curVertID = colResult.idVQuery;
		int32_t curMeshID = colResult.idTMQuery;

		int32_t collidedMeshID = collidingPt.intersectedMeshId;
		int32_t closestFaceId = collidingPt.closestSurfaceFaceId;
		TetMeshFEM* pCurTM = tMeshPtrs[curMeshID].get();
		TetMeshFEM* pIntersectedTM = tMeshPtrs[collidedMeshID].get();

		Vec3 faceNormal;
		pIntersectedTM->computeFaceNormal(closestFaceId, faceNormal);
		if (pointType == GAIA::ClosestPointOnTriangleType::AtInterior)
		{
			normal = faceNormal;
		}
		// all other cases
		else if (pointType != GAIA::ClosestPointOnTriangleType::NotFound)
		{
			Vec3 closestP_to_p;
			closestP_to_p = collidingPt.closestSurfacePt - pCurTM->vertex(curVertID);
			normal  = closestP_to_p.normalized();

			// the normal should point to the outside of the intersected mesh
			if (normal.dot(faceNormal) < 0)
			{
				normal = -normal;;
			}
		}

		/*switch (pointType)
		{
		case GAIA::ClosestPointOnTriangleType::AtA:
			surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(0, closestFaceId);
			surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(0, closestFaceId);
			pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);

			break;
		case GAIA::ClosestPointOnTriangleType::AtB:
			surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(1, closestFaceId);
			surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(1, closestFaceId);
			pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::AtC:
			surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(2, closestFaceId);
			surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(2, closestFaceId);
			pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::AtAB:
			pIntersectedTM->computeEdgeNormal(closestFaceId, 0, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::AtBC:
			pIntersectedTM->computeEdgeNormal(closestFaceId, 1, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::AtAC:
			pIntersectedTM->computeEdgeNormal(closestFaceId, 2, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::AtInterior:
			pIntersectedTM->computeFaceNormal(closestFaceId, normal);
			break;
		case GAIA::ClosestPointOnTriangleType::NotFound:
			return;
			break;
		default:
			return;
			break;
		}*/
	}

	inline void computeVFContactNormalTriMesh(IdType meshIdVertexSide, IdType vertexId, IdType meshIdFaceSide, IdType faceId, const Vec3& closestPt,
		std::vector<std::shared_ptr<TriMeshFEM>>& meshPtrs, ClosestPointOnTriangleType pointType, Vec3& normal)
	{
		const TriMeshFEM* pMeshVertexSide = meshPtrs[meshIdVertexSide].get();
		const TriMeshFEM* pMeshFaceSide = meshPtrs[meshIdFaceSide].get();

		Vec3 faceNormal = pMeshFaceSide->computeNormal(faceId);
		
		Vec3 closestP_to_p = pMeshVertexSide->vertex(vertexId) - closestPt;

		if (pointType == GAIA::ClosestPointOnTriangleType::AtInterior)
		{
			normal = faceNormal;

			// the contact normal should always points to the contact point, because it's not supposed to penetrate the mesh
			// if faceNormal.dot(closestP_to_) < 0 means the vertex is contacting the back face
			if (faceNormal.dot(closestP_to_p) < 0)
			{
				normal = -normal;
			}
		}
		// all other cases
		else if (pointType != GAIA::ClosestPointOnTriangleType::NotFound)
		{
			normal = closestP_to_p.normalized();
		}
	}

	inline void computeVFContactNormalTriMesh(const embree::Vec3fa& a, const embree::Vec3fa& b, const embree::Vec3fa& c, 
		const embree::Vec3fa& vertex, const embree::Vec3fa& closestPt, ClosestPointOnTriangleType pointType, Vec3& normal)
	{
		const embree::Vec3fa closestP_to_p = vertex - closestPt;
		const embree::Vec3fa faceNormal = embree::normalize_safe(embree::cross(b-a, c-a));

		if (pointType == GAIA::ClosestPointOnTriangleType::AtInterior)
		{
			normal << faceNormal.x, faceNormal.y, faceNormal.z;

			// the contact normal should always points to the contact point, because it's not supposed to penetrate the mesh
			// if faceNormal.dot(closestP_to_) < 0 means the vertex is contacting the back face
			if (embree::dot(faceNormal, closestP_to_p) < 0)
			{
				normal = -normal;
			}
		}
		// all other cases
		else if (pointType != GAIA::ClosestPointOnTriangleType::NotFound)
		{
			const embree::Vec3fa closestP_to_p_n = embree::normalize_safe(closestP_to_p);
			normal << closestP_to_p_n.x , closestP_to_p_n.y, closestP_to_p_n.z;
		}
	}
}

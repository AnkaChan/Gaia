#pragma once
#include "../common/math/vec2.h"
#include "../common/math/vec3.h"
#include "CollisionDetertionParameters.h"

#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"

#include <algorithm>

#define d_of(m, n, o, p) ((m.x - n.x) * (o.x - p.x) + (m.y - n.y) * (o.y - p.y) + (m.z - n.z) * (o.z - p.z))


namespace GAIA {

	struct TriMeshForCollision
	{
		TriMeshForCollision(TetMeshFEM* PTetMesh)
			: posBuffer(PTetMesh->vertices().data())
			, indexBuffer(PTetMesh->surfaceFacesTetMeshVIds().data())
			, numVertices(PTetMesh->numVertices())
			, numFaces(PTetMesh->numSurfaceFaces())
		{

		}
		TriMeshForCollision(TriMeshFEM* pTriMesh)
			: posBuffer(pTriMesh->positions().data())
			, indexBuffer(pTriMesh->facePos.data())
			, numVertices(pTriMesh->numVertices())
			, numFaces(pTriMesh->numFaces())
		{

		}

		IdType* getFaceVIds(IdType faceId) const
		{
			return (IdType*)(indexBuffer + 3 * faceId);
		}

		FloatingType* getVertex(IdType vertexId) const
		{
			return (FloatingType*)(posBuffer + 3 * vertexId);
		}

		Vec3 getVertexVec3(IdType vertexId) const
		{
			return Vec3(posBuffer[3 * vertexId], posBuffer[3 * vertexId + 1], posBuffer[3 * vertexId + 2]);
		}

		FloatingType* posBuffer;
		const IdType* indexBuffer;
		size_t numVertices;
		size_t numFaces;
	};

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

	// Copy from Godot Engine
	inline void get_closest_points_between_segments(const embree::Vec3fa& p_p0, const embree::Vec3fa& p_p1, const embree::Vec3fa& p_q0,
		const embree::Vec3fa& p_q1, embree::Vec3fa& r_ps, embree::Vec3fa& r_qt, FloatingType& s, FloatingType& t) {
		// Based on David Eberly's Computation of Distance Between Line Segments algorithm.

		embree::Vec3fa p = p_p1 - p_p0;
		embree::Vec3fa q = p_q1 - p_q0;
		embree::Vec3fa r = p_p0 - p_q0;

		FloatingType a = embree::dot(p, p);
		FloatingType b = embree::dot(p, q);
		FloatingType c = embree::dot(q, q);
		FloatingType d = embree::dot(p, r);
		FloatingType e = embree::dot(q, r);

		s = 0.0f;
		t = 0.0f;

		FloatingType det = a * c - b * b;
		if (det > CMP_EPSILON) {
			// Non-parallel segments
			FloatingType bte = b * e;
			FloatingType ctd = c * d;

			if (bte <= ctd) {
				// s <= 0.0f
				if (e <= 0.0f) {
					// t <= 0.0f
					s = (-d >= a ? 1 : (-d > 0.0f ? -d / a : 0.0f));
					t = 0.0f;
				}
				else if (e < c) {
					// 0.0f < t < 1
					s = 0.0f;
					t = e / c;
				}
				else {
					// t >= 1
					s = (b - d >= a ? 1 : (b - d > 0.0f ? (b - d) / a : 0.0f));
					t = 1;
				}
			}
			else {
				// s > 0.0f
				s = bte - ctd;
				if (s >= det) {
					// s >= 1
					if (b + e <= 0.0f) {
						// t <= 0.0f
						s = (-d <= 0.0f ? 0.0f : (-d < a ? -d / a : 1));
						t = 0.0f;
					}
					else if (b + e < c) {
						// 0.0f < t < 1
						s = 1;
						t = (b + e) / c;
					}
					else {
						// t >= 1
						s = (b - d <= 0.0f ? 0.0f : (b - d < a ? (b - d) / a : 1));
						t = 1;
					}
				}
				else {
					// 0.0f < s < 1
					FloatingType ate = a * e;
					FloatingType btd = b * d;

					if (ate <= btd) {
						// t <= 0.0f
						s = (-d <= 0.0f ? 0.0f : (-d >= a ? 1 : -d / a));
						t = 0.0f;
					}
					else {
						// t > 0.0f
						t = ate - btd;
						if (t >= det) {
							// t >= 1
							s = (b - d <= 0.0f ? 0.0f : (b - d >= a ? 1 : (b - d) / a));
							t = 1;
						}
						else {
							// 0.0f < t < 1
							s /= det;
							t /= det;
						}
					}
				}
			}
		}
		else {
			// Parallel segments
			if (e <= 0.0f) {
				s = (-d <= 0.0f ? 0.0f : (-d >= a ? 1 : -d / a));
				t = 0.0f;
			}
			else if (e >= c) {
				s = (b - d <= 0.0f ? 0.0f : (b - d >= a ? 1 : (b - d) / a));
				t = 1;
			}
			else {
				s = 0.0f;
				t = e / c;
			}
		}

		r_ps = (1 - s) * p_p0 + s * p_p1;
		r_qt = (1 - t) * p_q0 + t * p_q1;
	}

	inline FloatingType get_closest_distance_between_segments(const embree::Vec3fa& p_p0, const embree::Vec3fa& p_p1, const embree::Vec3fa& p_q0, const embree::Vec3fa& p_q1) {
		embree::Vec3fa ps;
		embree::Vec3fa qt;
		FloatingType s, t;
		get_closest_points_between_segments(p_p0, p_p1, p_q0, p_q1, ps, qt, s, t);
		embree::Vec3fa st = qt - ps;
		return embree::length(st);
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
			normal = closestP_to_p.normalized();

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
		const embree::Vec3fa faceNormal = embree::normalize_safe(embree::cross(b - a, c - a));

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
			normal << closestP_to_p_n.x, closestP_to_p_n.y, closestP_to_p_n.z;
		}
	}
}

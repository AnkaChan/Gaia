#pragma once

#include "../CyCodeBase/cyCore.h"
#include "../CyCodeBase/cyTriMesh.h"
#include "../CyCodeBase/cyMatrix.h"
#include "../Types/Types.h"
#define EDGE_UV_EPSILON 1e-6f



namespace GAIA {
	class TriData
	{
	public:
		Vec3 vs[3];	// triangle vertices
		Vec3 n;		// normal (not normalized)
		float d;		// plane equation scalar
		float dv[3];	// vertex distances to the other plane
		int outVerts;	// which vertices are on which side (0:below/inside or 1:above/outside), corresponds to the first 3 bits of the data
		int ix[4];		// vertex indices for the edge intersections. The vertex index on the opposite side is placed at [0] and [3]. 
						// The others are [1] and [2]. Possible values are {0,1,2,0}, {1,2,0,1}, and {2,0,1,2}
		float t[2];		// edge intersection position along the intersection line
		float f[2];		// edge intersection position along the edge
		int order;		// the order of edge intersections (0 or 1); if order = 0, the the tt is value for edge: (ix[0], ix[1]), (ix[2], ix[3]), respectively; 
		float tt[2];	// ordered t values
		Vec2 uv[2];	    // uv coordinates of the edge intersections
		Vec3 intersectionPoints[2];  // 3D position of the edge intersections

		void edgePosition(float f, int edgeId, Vec3& p) {
			int eV1 = edgeId;
			int eV2 = (edgeId + 1) % 3;

			p = vs[eV1] * (1 - f) + f * vs[eV2];
		}

		void baryPosition(Vec2 bary, Vec3& p) {

			p = vs[0] * bary[0] + vs[1] * bary[1] + vs[2] * (1- bary[0]- bary[1]);
		}
		
		TriData(const IdType* fVIds, const TVerticesMat& verts) {
			vs[0] = verts.col(fVIds[0]);
			vs[1] = verts.col(fVIds[1]);
			vs[2] = verts.col(fVIds[2]);
		}
		const Vec3& Normal() { return n; }
		int OutVerts() { return outVerts; }
		void SetPlane()
		{
			Vec3 e1 = vs[1] - vs[0];
			Vec3 e2 = vs[2] - vs[0];
			n = e1.cross(e2);
			d = -(n.dot(vs[0]));
		}
		// Returns true of there are no intersections
		bool ComputeDistances(TriData const& other)
		{
			for (int i = 0; i < 3; ++i) dv[i] = other.n.dot(vs[i]) + other.d;
			outVerts = (dv[0] >= 0) + ((dv[1] >= 0) << 1) + ((dv[2] >= 0) << 2);
			return (outVerts == 0) || (outVerts == 7);	// no intersection
		}
		void Barycentric(Vec3& p, float& u, float& v, float& w)
		{
			Vec3& a = vs[0];
			Vec3& b = vs[1]; 
			Vec3& c = vs[2];
			Vec3 v0 = b - a, v1 = c - a, v2 = p - a;
			float d00 = v0.dot(v0);
			float d01 = v0.dot(v1);
			float d11 = v1.dot(v1);
			float d20 = v2.dot(v0);
			float d21 = v2.dot(v1);
			float denom = d00 * d11 - d01 * d01;
			v = (d11 * d20 - d01 * d21) / denom;
			w = (d00 * d21 - d01 * d20) / denom;
			u = 1.0f - v - w;
		}
		// Compute intervals along the intersection line using the given dimension
		void ComputeIntervals(int dim)
		{
			const int s = ((outVerts & 4) ? ~outVerts : outVerts) & 3;

			ix[0] = s - 1;
			ix[1] = s % 3;
			ix[2] = (s + 1) % 3;
			ix[3] = ix[0];

			for (int i = 0, j = 0; i < 2; ++i, j += 2) {
				float p0 = vs[ix[j]][dim];
				float p1 = vs[ix[j + 1]][dim];
				float d0 = dv[ix[j]];
				float d1 = dv[ix[j + 1]];
				f[i] = d0 / (d0 - d1);
				t[i] = p0 + (p1 - p0) * f[i];
			}

			order = t[0] > t[1];
			tt[0] = t[order];
			tt[1] = t[1 - order];
		}
		float IntervalStart() const { return tt[0]; }
		float IntervalEnd() const { return tt[1]; }
		bool AreIntervalsSeparate(TriData const& b) const { return tt[1] < b.tt[0] || tt[0] > b.tt[1]; }
		void ComputeIntersections()
		{
			// const float f0[] = { 0.0f, f[0], 1 - f[0], 0.0f };
			// const float f1[] = { 0.0f, f[1], 1 - f[1], 0.0f };
			/*uv[0] << f0[ix[0] + 1], f0[ix[0]];
			uv[1] << f1[ix[2] + 1], f1[ix[2]];*/

			intersectionPoints[0] = vs[ix[0]] * (1.f-f[0]) + vs[ix[1]] * (f[0]);
			intersectionPoints[1] = vs[ix[2]] * (1.f-f[1]) + vs[ix[3]] * (f[1]);

		}
		Vec2 TriPos(float _t) const { return uv[0] + (uv[1] - uv[0]) * ((_t - t[0]) / (t[1] - t[0])); }

	};

	class alignas(32) TriTriIntersection
	{
	public:
		Vec2              uv[2][2];	   // of the two ends, iFace x iIntersection
		float             fs[2];       // of the two ends, edge intersection position along the edge. only one parameter is recored, another can be obtained by substract by 1
		int meshIds[2];
		struct IDs {
			uint8_t s : 1;	// side index of the end
			uint8_t e : 2;	// edge index of the end
			uint8_t downIntersection : 1;
			uint8_t unused : 4;
		};
		IdType            fid[2];	// of the two faces
		uint8_t           outV[2];	// 3 bits per face, showing which vertices are below the plane of the other (1:above/outside, 0:below/inside).
		IDs               ix[2];	// End of the intersection line segment, from left (smaller t) to right (larger t)
		// It shows which edge marks the two ends of the intersection line segment. They can be any one of the 3 edges of the two triangles (6 in total).
		// s indicates which triangle (i.e. 0 or 1) and e is the edge index of that triangle, such that edge 0 is from vertex 0 to 1, edge 1 is from vertex 1 to 2, 
		// and edge 2 is from vertex 2 to 0.

		float t[2];	// position along the intersecting edge of the two ends

		void setMeshIds(IdType mesh1, IdType mesh2) {
			meshIds[0] = mesh1;
			meshIds[1] = mesh2;
		}

		inline Vec3 get3DPosition(int iFace, int iIntersect, IdType* fVIds, const TVerticesMat& verts ) const {
			Vec3 pt = verts.col(fVIds[0]) * uv[iFace][iIntersect](0) + verts.col(fVIds[1]) * uv[iFace][iIntersect](1)
				+ verts.col(fVIds[2]) * (1.f - uv[iFace][iIntersect](0)- uv[iFace][iIntersect](1));

			return pt;
		}

		int getIFace(int meshId, int faceId) const {
			int iFace = -1;
			if (faceId == fid[0] && meshId == meshIds[0])
				iFace = 0;
			else if (faceId == fid[1] && meshId == meshIds[1])
				iFace = 1;

			return iFace;
		}

		// twoEndsMask: first bit : whether end point 1 is and edge intersection of this face
		//              second bit: whether end point 2 is and edge intersection of this face
		bool isStartingEdge(int iFace, int & twoEndsMask) const {
			twoEndsMask = 0;
			if (ix[0].s == iFace)
			{
				twoEndsMask = twoEndsMask ^ 1;
			}
			if (ix[1].s == iFace)
			{
				twoEndsMask = twoEndsMask ^ 2;
			}

			return twoEndsMask;
		}

		bool Intersect(IdType fVId1, IdType fVId2, const IdType* fVIds1, const IdType* fVIds2, const TVerticesMat & verts1, const TVerticesMat& verts2);
		IdType maxIndex(const Vec3& v) {
			int maxId = v[0] > v[1] ? 0 : 1;
			maxId = v[maxId] > v[2] ? maxId : 2;
			return maxId;
		}
	};
	inline bool GAIA::TriTriIntersection::Intersect(IdType fVId1, IdType fVId2, const IdType* fVIds1, const IdType* fVIds2, const TVerticesMat& verts1, const TVerticesMat& verts2)
	{
		for (int iFV1 = 0; iFV1 < 3; iFV1++)
		{
			for (int iFV2 = 0; iFV2 < 3; iFV2++) {
				if (fVIds1[iFV1] == fVIds2[iFV2])
				{
					return false;
				}
			}
		}


		TriData ta(fVIds1, verts1);
		TriData tb(fVIds2, verts2);
		TriData* tris[2] = {
			&ta, &tb
		};

		ta.SetPlane();
		if (tb.ComputeDistances(ta)) return false;

		tb.SetPlane();
		if (ta.ComputeDistances(tb)) return false;

		const Vec3 d = ta.Normal().cross(tb.Normal());	// intersection line direction (not normalized)
		const Vec3 d_abs = d.cwiseAbs();
		const int dim = maxIndex(d_abs);
		const float len2 = ta.Normal().squaredNorm() * tb.Normal().squaredNorm();
		const float d_max2 = d_abs[dim] * d_abs[dim];
		if (d_max2 < 1e-5 * len2) return false;	// the triangles are parallel.

		ta.ComputeIntervals(dim);
		tb.ComputeIntervals(dim);

		if (ta.AreIntervalsSeparate(tb)) return false; // the intervals don't overlap

		// There is an intersection. We must now compute what it is

		ta.ComputeIntersections();
		tb.ComputeIntersections();

		fid[0] = fVId1;
		fid[1] = fVId2;

		// whichever has the larger IntervalStart, makes the left(less) side of the line seg
		int i0 = ta.IntervalStart() < tb.IntervalStart();
		// i0 == 0, means, ta makes the left end
		ix[0].s = i0;
		// order = 0, means the tt is value for edge: (ix[0], ix[1]), (ix[2], ix[3])
		ix[0].e = tris[i0]->ix[(tris[i0]->order) * 2];
		Vec3& p0 = tris[i0]->intersectionPoints[tris[i0]->order];
		float w;
		ta.Barycentric(p0, uv[0][0](0), uv[0][0](1), w);
		tb.Barycentric(p0, uv[1][0](0), uv[1][0](1), w);
		fs[0] = tris[i0]->f[(tris[i0]->order)];
		// compute whether it's an up or down intersection
		int targetVert = tris[i0]->ix[(tris[i0]->order) * 2 + 1];
		ix[0].downIntersection = tris[i0]->outVerts & (1 << targetVert) ? 1 :0;

		int i1 = ta.IntervalEnd() > tb.IntervalEnd();
		ix[1].s = i1;
		// order = 0, means the tt is value for edge: (ix[0], ix[1]), (ix[2], ix[3])
		ix[1].e = tris[i1]->ix[(tris[i1]->order^1) * 2];
		Vec3& p1 = tris[i1]->intersectionPoints[tris[i1]->order ^ 1];
		ta.Barycentric(p1, uv[0][1](0), uv[0][1](1), w);
		tb.Barycentric(p1, uv[1][1](0), uv[1][1](1), w);
		fs[1] = tris[i1]->f[(tris[i1]->order ^ 1)];
		targetVert = tris[i1]->ix[(tris[i1]->order ^ 1) * 2 + 1];
		ix[1].downIntersection = tris[i1]->outVerts & (1 << targetVert) ? 1 : 0;

		outV[0] = ta.OutVerts();
		outV[1] = tb.OutVerts();
		// ta.tt[0] corrsponds to edge: (ix[0], ix[1]) if order == 0llll
		// ta.tt[0] corrsponds to edge: (ix[2], ix[3]) if order == 1
		// (ix[order *2], ix[order*2 + 1])
		return true;
	}

}


//namespace cy {
//	typedef int FaceID;
//	class alignas(32) TriTriIntersection
//	{
//	public:
//		struct IDs {
//			uint8_t s : 1;	// side index of the end
//			uint8_t e : 2;	// edge index of the end
//			uint8_t unused : 5;
//		};
//		Vec2f             uv[2];	// of the two ends
//		FaceID            fid[2];	// of the two sides
//		uint8_t           outV[2];	// 3 bits per face, showing which vertices are below the plane of the other (1:below, 0:above).
//		IDs               ix[2];	// of the two ends
//
//		// It shows which edge marks the two ends of the intersection line segment. They can be any one of the 3 edges of the two triangles (6 in total).
//		// s indicates which triangle (i.e. 0 or 1) and e is the edge index of that triangle, such that edge 0 is from vertex 0 to 1, edge 1 is from vertex 1 to 2, 
//		// and edge 2 is from vertex 2 to 0.
//		float t[2];	// position along the intersecting edge of the two ends
//		/*template <typename MeshTypeA, typename MeshTypeB>
//		bool Intersect(MeshTypeA const& ma, FaceID a, MeshTypeB const& mb, FaceID b)
//		{
//			Vec3f va[3], vb[3];
//			ma.GetVerts(va, a);
//			mb.GetVerts(vb, b);
//			return Intersect(a, b, va, vb);
//		}
//
//		template <typename MeshType>
//		bool Intersect(MeshType const& m, FaceID a, FaceID b)
//		{
//			if (a > b) Swap(a, b);
//			return Intersect(m, a, m, b);
//		}*/
//
//		bool operator < (TriTriIntersection const& i) const { int j = fid[0] == i.fid[0]; return fid[j] < i.fid[j]; }
//		FaceID  Face(int side) const { return fid[side]; }
//		int     FaceSide(FaceID f) const { return fid[1] == f; }
//		uint8_t OutsideVerts(int side) const { return outV[side]; }
//		int     IntervalSide(int end)  const { return ix[end].s; }
//		FaceID  IntervalFace(int end)  const { return fid[ix[end].s]; }
//		int     IntervalEdgeIndex(int end)  const { return ix[end].e; }
//		float   IntervalParam(int end)  const { return t[end]; }
//		Vec2f   IntervalUV(int end)  const { return UV(t[end], ix[end].e); }
//		bool    IntervalIsFrontHit(int end)  const { return (outV[ix[end].s] & (1 << ix[end].e)) > 0; }	// Is the edge intersection of the interval end going from outside to inside, using edge vertex orders 0->1, 1->2, or 2->0.
//		int     IntervalInVertIndex(int end)  const { int vix = IntervalEdgeIndex(end); return IntervalIsFrontHit(end) ? (vix + 1) % 3 : vix; }	// returns the index of the face vertex that is inside
//		FaceID  IntervalOtherFace(int end)  const { return fid[1 - ix[end].s]; }
//		Vec2f   IntervalOtherUV(int end)  const { return uv[end]; }
//		bool    FindIntervalEnd(int& end, FaceID f, int edgeIx) const
//		{
//			for (end = 0; end < 2; ++end) {
//				if (fid[ix[end].s] == f && ix[end].e == edgeIx) return true;
//			}
//			// In case of edge-edge intersection, we might be looking at the wrong edge.
//			for (end = 0; end < 2; ++end) {
//				if (fid[1 - ix[end].s] == f) {
//					float dif[3] = {
//						uv[end].y,
//						1 - uv[end].Sum(),
//						uv[end].x,
//					};
//					if (std::abs(dif[edgeIx]) < EDGE_UV_EPSILON) return true;
//				}
//			}
//			return false;
//		}
//		bool FindIntervalSideEnd(int& end, int side, int edgeIx) const
//		{
//			for (end = 0; end < 2; ++end) {
//				if (ix[end].s == side && ix[end].e == edgeIx) return true;
//			}
//			// In case of edge-edge intersection, we might be looking at the wrong edge.
//			for (end = 0; end < 2; ++end) {
//				if (ix[end].s != side) {
//					float dif[3] = {
//						uv[end].y,
//						1 - uv[end].Sum(),
//						uv[end].x,
//					};
//					if (std::abs(dif[edgeIx]) < EDGE_UV_EPSILON) return true;
//				}
//			}
//			return false;
//		}
//		bool FindIntervalEnd(int& end, FaceID f, int edgeIx, int side, bool selfCollision) const
//		{
//			return selfCollision ? FindIntervalEnd(end, f, edgeIx) : FindIntervalSideEnd(end, side, edgeIx);
//		}
//		//bool IsVertInside( FaceID f, int vertIx ) const { return ( sides[ fid[1]==f ] >> vertIx ) & 1; }
//		bool IsVertInside(int side, int vertIx) const { assert(side == 0 || side == 1); return (outV[side] >> vertIx) & 1; }
//		void Invalidate() { ix[0].e = 3; }
//		bool IsValid() const { return ix[0].e != 3; }
//
//		static Vec2f UV(float t, int edgeIx) { const float f[] = { 0,t,1 - t,0 }; return Vec2f(f[edgeIx + 1], f[edgeIx]); }
//
//		//Vec3f IntervalPoint(TMeshBase const& m, int end) const
//		//{
//		//	if (ix[end].s == 0) {
//		//		return IntervalPointB(m, end);
//		//	}
//		//	else {
//		//		return IntervalPointA(m, end);
//		//	}
//		//}
//
//		//Vec3f IntervalPointA(TMeshBase const& m, int end) const
//		//{
//		//	FaceID fi = fid[1 - ix[end].s];
//		//	Vec3f v[3];
//		//	m.GetVerts(v, fi);
//		//	Vec2f u = uv[end];
//		//	return v[0] * (1 - u.x - u.y) + v[1] * u.x + v[2] * u.y;
//		//}
//		//Vec3f IntervalPointB(TMeshBase const& m, int end) const
//		//{
//		//	FaceID fi = fid[ix[end].s];
//		//	Vec3f v[3];
//		//	m.GetVerts(v, fi);
//		//	float tv = float(t[end]);
//		//	Vec3f v0 = v[ix[end].e];
//		//	Vec3f v1 = v[(ix[end].e + 1) % 3];
//		//	return v0 + (v1 - v0) * tv;
//		//}
//		/*
//		void Print()
//		{
//			printf("%d %d --> [ %d %d %f ] [ %d %d %f ]\n", (int)fid[0], (int)fid[1], pt[0].f_ix, pt[0].e_ix, float(t[0]), pt[1].f_ix, pt[1].e_ix, float(t[1]) );
//		}
//		*/
//		bool Intersect(FaceID a, FaceID b, Vec3f const va[3], Vec3f const vb[3]);
//		static float EdgePos(float& t, float& f, float p0, float p1, float d0, float d1)
//		{
//			f = d0 / (d1 - d0);
//			t = p0 + (p1 - p0) * f;
//		}
//	};
//
//	//-------------------------------------------------------------------------------
//
//	inline bool TriTriIntersection::Intersect(FaceID a, FaceID b, Vec3f const va[3], Vec3f const vb[3])
//	{
//		class TriData
//		{
//			Vec2f uv[2];	// uv coordinates of the edge intersections
//			Vec3f const* v;	// triangle vertices
//			Vec3f n;		// normal (not normalized)
//			float d;		// plane equation scalar
//			float dv[3];	// vertex distances to the other plane
//			int outVerts;	// which vertices are on which side (0:below/inside or 1:above/outside), corresponds to the first 3 bits of the data
//			int ix[4];		// vertex indices for the edge intersections. The vertex index on the opposite side is placed at [0] and [3]. 
//							// The others are [1] and [2]. Possible values are {0,1,2,0}, {1,2,0,1}, and {2,0,1,2}
//			float t[2];		// edge intersection position along the intersection line
//			float f[2];		// edge intersection position along the edge
//			int order;		// the order of edge intersections (0 or 1)
//			float tt[2];	// ordered t values
//		public:
//			void Init(Vec3f const* _v) { v = _v; }
//			Vec3f const& Normal() const { return n; }
//			int OutVerts() const { return outVerts; }
//			void SetPlane()
//			{
//				const Vec3f e1 = v[1] - v[0];
//				const Vec3f e2 = v[2] - v[0];
//				n = e1 ^ e2;
//				d = -(n % v[0]);
//			}
//			// Returns true of there are no intersections
//			bool ComputeDistances(TriData const& other)
//			{
//				for (int i = 0; i < 3; ++i) dv[i] = other.n % v[i] + other.d;
//				outVerts = (dv[0] >= 0) + ((dv[1] >= 0) << 1) + ((dv[2] >= 0) << 2);
//				return (outVerts == 0) || (outVerts == 7);	// no intersection
//			}
//			// Compute intervals along the intersection line using the given dimension
//			void ComputeIntervals(int dim)
//			{
//				const int s = ((outVerts & 4) ? ~outVerts : outVerts) & 3;
//
//				ix[0] = s - 1;
//				ix[1] = s % 3;
//				ix[2] = (s + 1) % 3;
//				ix[3] = ix[0];
//
//				for (int i = 0, j = 0; i < 2; ++i, j += 2) {
//					float p0 = v[ix[j]][dim];
//					float p1 = v[ix[j + 1]][dim];
//					float d0 = dv[ix[j]];
//					float d1 = dv[ix[j + 1]];
//					f[i] = d0 / (d0 - d1);
//					t[i] = p0 + (p1 - p0) * f[i];
//				}
//
//				order = t[0] > t[1];
//				tt[0] = t[order];
//				tt[1] = t[1 - order];
//			}
//			float IntervalStart() const { return tt[0]; }
//			float IntervalEnd() const { return tt[1]; }
//			bool AreIntervalsSeparate(TriData const& b) const { return tt[1] < b.tt[0] || tt[0] > b.tt[1]; }
//			void ComputeIntersections()
//			{
//				const float f0[] = { 0.0f, f[0], 1 - f[0], 0.0f };
//				const float f1[] = { 0.0f, f[1], 1 - f[1], 0.0f };
//				uv[0].Set(f0[ix[0] + 1], f0[ix[0]]);
//				uv[1].Set(f1[ix[2] + 1], f1[ix[2]]);
//			}
//
//		};
//
//		TriData tdata[2];
//		TriData& ta = tdata[0];
//		TriData& tb = tdata[1];
//		ta.Init(va);
//		tb.Init(vb);
//
//		ta.SetPlane();
//		if (tb.ComputeDistances(ta)) return false;
//
//		tb.SetPlane();
//		if (ta.ComputeDistances(tb)) return false;
//
//		const Vec3f d = ta.Normal() ^ tb.Normal();	// intersection line direction (not normalized)
//		const Vec3f d_abs = d.Abs();
//		const int dim = d_abs.MaxIndex();
//		const float len2 = ta.Normal().LengthSquared() * tb.Normal().LengthSquared();
//		const float d_max2 = d_abs[dim] * d_abs[dim];
//		if (d_max2 < 1e-5 * len2) return false;	// the triangles are parallel.
//
//		ta.ComputeIntervals(dim);
//		tb.ComputeIntervals(dim);
//
//		if (ta.AreIntervalsSeparate(tb)) return false; // the intervals don't overlap
//
//		// There is an intersection. We must now compute what it is
//
//		ta.ComputeIntersections();
//		tb.ComputeIntersections();
//
//		fid[0] = a;
//		fid[1] = b;
//
//		int i0 = ta.IntervalStart() < tb.IntervalStart();
//	}
//
//}
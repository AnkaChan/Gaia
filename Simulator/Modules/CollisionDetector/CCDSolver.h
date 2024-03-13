#ifndef __CCD_SOLVER__
#define __CCD_SOLVER__

#include "../CyCodeBase/cyCore.h"
#include "../CyCodeBase/cyVector.h"
#include "../CyCodeBase/cyPolynomial.h"
#include <cmath>

#if defined(__INTEL_COMPILER) || defined(__clang__) || defined(__GNUC__)
# define nodefault default: __builtin_unreachable()
#elif defined(_MSC_VER)
# define nodefault default: __assume(0)
#else
# define nodefault default: 
#endif



namespace cy {

	//-------------------------------------------------------------------------------

	template <typename T, typename S> inline T _MultSign(T v, S sign) { return v * (sign < 0 ? T(-1) : T(1)); }	//!< Multiplies the given value with the given sign

	//-------------------------------------------------------------------------------

	//template <typename ftype>
	//inline int QuadraticRoots(ftype roots[2], ftype const poly[3], ftype x0, ftype x1)
	//{
	//	const ftype c = poly[0];
	//	const ftype b = poly[1];
	//	const ftype a = poly[2];
	//	const ftype delta = b * b - 4 * a * c;
	//	if (delta >= 0) {
	//		const ftype d = Sqrt(delta);
	//		const ftype q = ftype(-0.5) * (b + _MultSign(d, b));
	//		const ftype rv0 = q / a;
	//		const ftype rv1 = c / q;
	//		const ftype r0 = Min(rv0, rv1);
	//		const ftype r1 = Max(rv0, rv1);
	//		int r0i = (r0 >= x0) & (r0 <= x1);
	//		int r1i = (r1 >= x0) & (r1 <= x1);
	//		roots[0] = r0;
	//		roots[r0i] = r1;
	//		return r0i + r1i;
	//	}
	//	return 0;
	//}

	//-------------------------------------------------------------------------------

	template <typename ftype>
	inline int CubicRootsAnalytic(ftype roots[3], ftype const poly[4], ftype x0, ftype x1)
	{
		if (poly[3] != 0) {
			const ftype a = poly[2] / poly[3];
			const ftype b = poly[1] / poly[3];
			const ftype c = poly[0] / poly[3];
			const ftype q = a * a - 3 * b;
			const ftype r = a * a * a - ftype(4.5) * a * b + ftype(13.5) * c;
			const ftype q3 = q * q * q;
			const ftype r2 = r * r;
			const ftype r2_q3 = r2 - q3;
			if (r2_q3 >= 0) {
				// single real root
				const ftype e = _MultSign(cbrt(std::abs(r) + Sqrt(r2_q3)), -r);
				const ftype f = e != 0 ? q / e : 0;
				const ftype r0x3 = e + f - a;
				roots[0] = r0x3 / 3;
				return (r0x3 >= 3 * x0 && r0x3 <= 3 * x1);
			}
			else {
				// three real roots
				const ftype sq = Sqrt(q);
				const ftype sq2 = 2 * sq;
				if (-sq2 - a > 3 * x1) return 0;
				if (sq2 - a < 3 * x0) return 0;
				const ftype theta = std::acos(r / (q * sq));
				const ftype theta_3 = theta / 3;
				ftype _rx3[3];
				_rx3[0] = -sq2 * std::cos(theta_3) - a;
				_rx3[1] = -sq2 * std::cos(theta_3 + 2 * Pi<ftype>() / 3) - a;
				_rx3[2] = -sq2 * std::cos(theta_3 + 4 * Pi<ftype>() / 3) - a;
				ftype rx3[3];
				Sort3<true>(rx3, _rx3);
				int numRoots = 0;
				if (rx3[0] > x1 * 3) return 0;
				if (rx3[0] >= x0 * 3) { roots[0] = rx3[0] / 3; numRoots = 1; }
				if (rx3[1] > x1 * 3) return numRoots;
				if (rx3[1] >= x0 * 3) roots[numRoots++] = rx3[1] / 3;
				if (rx3[2] > x1 * 3) return numRoots;
				if (rx3[2] >= x0 * 3) roots[numRoots++] = rx3[2] / 3;
				return numRoots;
			}
		}
		else return QuadraticRoots<ftype>(roots, poly, x0, x1);
	}

	//-------------------------------------------------------------------------------

	// returns the number of t values in [0,1] when the points are planar
	template<typename DType>
	inline int MovingPointsPlanar(DType t[3], Vec3<DType> const x0[4], Vec3<DType> const x1[4])
	{
		Vec3<DType> dx0 = x1[0] - x0[0];
		Vec3<DType> dx1 = x1[1] - x0[1];
		Vec3<DType> dx2 = x1[2] - x0[2];
		Vec3<DType> dx3 = x1[3] - x0[3];

		Vec3<DType> e1 = x0[1] - x0[0];
		Vec3<DType> e2 = x0[2] - x0[0];
		Vec3<DType> e3 = x0[3] - x0[0];

		Vec3<DType> de0 = dx3 - dx0;
		Vec3<DType> de1 = dx1 - dx0;
		Vec3<DType> de2 = dx2 - dx0;

		Vec3<DType> e1xe2 = e1 ^ e2;
		Vec3<DType> de1xde2 = de1 ^ de2;
		Vec3<DType> dde1xe2 = (e1 ^ de2) + (de1 ^ e2);

		DType polynomial[4];
		polynomial[0] = e3 % e1xe2;
		polynomial[1] = (de0 % e1xe2) + (e3 % dde1xe2);
		polynomial[2] = (e3 % de1xde2) + (de0 % dde1xe2);
		polynomial[3] = de0 % de1xde2;

		//return CubicRootsAnalytic<DType>(t, polynomial, 0, 1);
		return CubicRoots<DType>(t, polynomial, 0, 1);
	}

	//-------------------------------------------------------------------------------
	template<typename DType>
	inline void InterpolatedPoints2D(Vec2<DType> y[4], DType t, Vec3<DType> const x0[4], Vec3<DType> const x1[4])
	{
		Vec3<DType> xx0 = x0[0] + (x1[0] - x0[0]) * t;
		Vec3<DType> xx1 = x0[1] + (x1[1] - x0[1]) * t;
		Vec3<DType> xx2 = x0[2] + (x1[2] - x0[2]) * t;
		Vec3<DType> xx3 = x0[3] + (x1[3] - x0[3]) * t;
		Vec3<DType> nrm = (xx1 - xx0) ^ (xx2 - xx0);
		int maxCmp = nrm.Abs().MaxComp();
		switch (maxCmp) {
		case 0: y[0] = xx0.YZ(); y[1] = xx1.YZ(); y[2] = xx2.YZ(); y[3] = xx3.YZ(); break;
		case 1: y[0] = xx0.XZ(); y[1] = xx1.XZ(); y[2] = xx2.XZ(); y[3] = xx3.XZ(); break;
		case 2: y[0] = xx0.XY(); y[1] = xx1.XY(); y[2] = xx2.XY(); y[3] = xx3.XY(); break;
		//default: break;
		nodefault;
		}
	}

	template<typename DType>
	inline void Barycentric(const Vec3<DType>& p, const  Vec3<DType>& a, const Vec3<DType>& b, const Vec3<DType>& c, Vec3<DType>& intersectionBarys)
	{
		Vec3<DType> v0 = b - a, v1 = c - a, v2 = p - a;
		DType d00 = v0.Dot(v0);
		DType d01 = v0.Dot(v1);
		DType d11 = v1.Dot(v1);
		DType d20 = v2.Dot(v0);
		DType d21 = v2.Dot(v1);
		DType denom = d00 * d11 - d01 * d01;
		intersectionBarys.y = (d11 * d20 - d01 * d21) / denom;
		intersectionBarys.z = (d00 * d21 - d01 * d20) / denom;
		intersectionBarys.x = 1.0f - intersectionBarys.y - intersectionBarys.z;
	}

	//-------------------------------------------------------------------------------

	// Moving triangle with vertices x and a moving point with position p
	template<typename DType>
	inline bool IntersectContinuousTriPoint(DType& tt, Vec3<DType> const x[2][3], Vec3<DType> const p[2], Vec3<DType> &intersectionBarys)
	{
		Vec3<DType> xx0[4] = { x[0][0], x[0][1], x[0][2], p[0] };
		Vec3<DType> xx1[4] = { x[1][0], x[1][1], x[1][2], p[1] };

		DType t[3];
		int numt = MovingPointsPlanar(t, xx0, xx1);

		for (int i = 0; i < numt; ++i) {
			// check if the plane intersection is inside the triangle
			Vec2<DType> y[4];
			InterpolatedPoints2D<DType>(y, t[i], xx0, xx1);
			y[0] -= y[3];
			y[1] -= y[3];
			y[2] -= y[3];
			bool y01p = std::signbit(y[0] ^ y[1]);
			bool y12p = std::signbit(y[1] ^ y[2]);
			bool y20p = std::signbit(y[2] ^ y[0]);
			if (y01p == y12p && y01p == y20p) {
				tt = t[i];

				Barycentric(p[0] + (p[1] - p[0]) * tt,
					x[0][0] + (x[1][0] - x[0][0]) * tt,
					x[0][1] + (x[1][1] - x[0][1]) * tt,
					x[0][2] + (x[1][2] - x[0][2]) * tt,
					intersectionBarys
				);
				return true;
			}
		}

		return false;
	}

	template<typename DType>
	inline bool IntersectContinuousEdgeEdge(DType& tt, Vec3<DType> const e0[2][2], Vec3<DType> const e1[2][2])
	{
		Vec3<DType> xx0[4] = { e0[0][0], e0[0][1], e1[0][0], e1[0][1] };
		Vec3<DType> xx1[4] = { e0[1][0], e0[1][1], e1[1][0], e1[1][1] };
		DType t[3];
		int numt = MovingPointsPlanar(t, xx0, xx1);

		for (int i = 0; i < numt; ++i) {
			// check if the edges overlap
			Vec2<DType> y[4];
			InterpolatedPoints2D(y, t[i], xx0, xx1);
			bool t012p = std::signbit((y[1] - y[0]) ^ (y[2] - y[0]));
			bool t013p = std::signbit((y[1] - y[0]) ^ (y[3] - y[0]));
			if (t012p == t013p) continue;
			bool t230p = std::signbit((y[3] - y[2]) ^ (y[0] - y[2]));
			bool t231p = std::signbit((y[3] - y[2]) ^ (y[1] - y[2]));
			if (t230p != t231p)
			{
				tt = t[i];

				return true;
			}
		}
		return false;
	}
};

#endif
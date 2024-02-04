#pragma once

#include "Eigen/core"
#include "Eigen/Dense"
#include "Eigen/Sparse"

#include <iostream>

#define GRAVITY_AXIS 1

#define POINT_VEC_DIMS 3
#define F_FLATTEN_DIMS 9
#define X_FLATTEN_DIMS 12

#define VEC_BLOCK_DIM3x1_AT(v, iV) (v.block<POINT_VEC_DIMS, 1>(POINT_VEC_DIMS * iV, 0))

#define SQR(x) ((x) * (x))
#define CUBE(x) ((x) * (x) * (x))

#define CMP_EPSILON 0.00001f
#define CMP_EPSILON2 (CMP_EPSILON * CMP_EPSILON)

#define CMP_NORMALIZE_TOLERANCE 0.000001f
#define CMP_POINT_IN_PLANE_EPSILON 0.00001f

namespace GAIA {
	typedef float FloatingType;
	typedef const float CFloatingType;
	//typedef double FloatingType;
	typedef int32_t IdType;

	// deformation gradient computation
	typedef Eigen::Matrix<FloatingType, 2, 1> Vec2;
	typedef Eigen::Matrix<IdType, 2, 1> Vec2I;
	typedef Eigen::Matrix<FloatingType, 3, 1> Vec3;
	typedef Eigen::Matrix<FloatingType, 4, 1> Vec4;
	typedef Eigen::Matrix<double, 3, 1> Vec3d;
	typedef Eigen::Matrix<IdType, 3, 1> Vec3I;
	typedef Eigen::Matrix<FloatingType, F_FLATTEN_DIMS, 1> Vec9;
	typedef Eigen::Matrix<FloatingType, X_FLATTEN_DIMS, 1> Vec12;
	typedef Eigen::Matrix<FloatingType, POINT_VEC_DIMS, POINT_VEC_DIMS> Mat3;
	typedef Eigen::Matrix<FloatingType, 9, 9> Mat9;
	typedef Eigen::Matrix<FloatingType, 12, 12> Mat12;
	typedef Eigen::Matrix<FloatingType, 4, 4> Mat4;
	typedef Eigen::Matrix<FloatingType, 2, 2> Mat2;
	typedef Eigen::Matrix<FloatingType, 3, 2> Mat3x2;
	typedef Eigen::Matrix<double, POINT_VEC_DIMS, POINT_VEC_DIMS> Mat3d;

	// tet mesh data structure
	typedef std::array<IdType, 4> TTetVIdsArr;
	typedef std::array<FloatingType, POINT_VEC_DIMS> TVertexCoordsArr;
	typedef std::array<FloatingType, 2> TVertexUVArr;
	typedef Eigen::Matrix<IdType, 4, Eigen::Dynamic> TTetIdsMat;
	typedef Eigen::Matrix<IdType, 2, Eigen::Dynamic> Mat2xI;
	typedef Eigen::Matrix<FloatingType, POINT_VEC_DIMS, Eigen::Dynamic> TVerticesMat;
	typedef Eigen::Matrix<FloatingType, 2, Eigen::Dynamic> TVerticesUVMat;
	typedef Eigen::Matrix<FloatingType, Eigen::Dynamic, 1> VecDynamic;
	typedef Eigen::Matrix<IdType, Eigen::Dynamic, 1> VecDynamicI;

	typedef Eigen::Matrix<FloatingType, 9, Eigen::Dynamic> Mat9x;
	typedef Eigen::Matrix<FloatingType, 12, Eigen::Dynamic> Mat12x;

	// on GPU the size of bool is 8 bits
	typedef Eigen::Matrix<int8_t, Eigen::Dynamic, 1> VecDynamicBool;
	typedef Eigen::Matrix<IdType, 3, Eigen::Dynamic> FaceVIdsMat;

	using Vec4BlockI = Eigen::Block<TTetIdsMat, 4, 1>;
	using Vec3Block = Eigen::Block<TVerticesMat, 3, 1>;
	using ConstVec3Block = Eigen::Block<const TVerticesMat, 3, 1>;

	// sparse matrix
	typedef Eigen::SparseMatrix<FloatingType> SpMat;
	typedef Eigen::Triplet<FloatingType> Triplet;

	namespace Utility {

		inline void printVecInfo(VecDynamic input, std::string name) {
			VecDynamic inputAbs = input.cwiseAbs();
			std::cout << name << ": " << "max " << inputAbs.maxCoeff() << " | min " << inputAbs.minCoeff() << " | mean " << inputAbs.mean() << std::endl;;
		}
		inline bool isZeroApprox(double s) {
			return abs(s) < CMP_EPSILON;
		}
	}



}
#pragma once

#include "../Types/Types.h"
#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"

#include "../CollisionDetector/ClothContactDetector.h"

#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif

namespace GAIA {
	//For newton solverDirect
#ifdef USE_DOUBLE
	using NFloatingType = double;
#else
	using NFloatingType = float;
#endif
	using NCFloatingType = const NFloatingType;
	using NSpMat = Eigen::SparseMatrix<NFloatingType>;
	using NVecDynamic = Eigen::Matrix<NFloatingType, Eigen::Dynamic, 1>;
	using NMat12 = Eigen::Matrix<NFloatingType, 12, 12>;
	using NMat9 = Eigen::Matrix<NFloatingType, 9, 9>;
	using NMat3 = Eigen::Matrix<NFloatingType, 3, 3>;
	using NVec9 = Eigen::Vector<NFloatingType, 9>;
	using NVec12 = Eigen::Vector<NFloatingType, 12>;
	using NVec3 = Eigen::Vector<NFloatingType, 3>;
	using NTriplet = Eigen::Triplet<NFloatingType>;
	using NTVerticesMat = Eigen::Matrix<NFloatingType, POINT_VEC_DIMS, Eigen::Dynamic>;
	using NMat = Eigen::Matrix<NFloatingType, Eigen::Dynamic, Eigen::Dynamic>;
	using NMati = Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic>;

#ifdef USE_MKL
	using DirectSolverType = Eigen::PardisoLDLT<NSpMat>;
	// using DirectSolverType = Eigen::PardisoLU<NSpMat>;
	using CGSolverType = Eigen::ConjugateGradient<NSpMat, Eigen::Lower | Eigen::Upper>;
#else
	using DirectSolverType = Eigen::SimplicialLDLT<NSpMat>;
	// using DirectSolverType = Eigen::SparseLU<NSpMat>;
	using CGSolverType = Eigen::ConjugateGradient<NSpMat, Eigen::Lower | Eigen::Upper>;
#endif

	struct BaseNewtonAssembler {

		size_t numAllVertices{};
		size_t numAllEdges{};

		std::vector<NTriplet> newtonHessianTripletsElasticity{};
		std::vector<NTriplet> newtonHessianTripletsBending{};

		// pointers to the diagonal blocks of the Hessian matrix
		std::vector<std::vector<NFloatingType*>> diagonalHessianBlockPtrs{};
		// pointers to the off-diagonal blocks of the Hessian matrix, not necessarily corresponding to the edges
		std::vector<std::vector<NFloatingType*>> offDiagonalHessianBlockPtrs{};

		// pointers to the Hessian blocks corresponds to collision, not necessarily corresponding to the edges
		// each VF collisoin has 4*3*3=36 elements
		std::vector<std::vector<NFloatingType*>> VFCollisionHessianBlockPtrs{};
		// each VF collisoin has 2*3*3=18 elements, it only corresponds to the edge of this side
		// the vertices from the other side are handled by the other side
		std::vector<std::vector<NFloatingType*>> EECollisionHessianBlockPtrs{};

		std::vector<NFloatingType> bendingHessianValues{};


		NFloatingType newtonEnergy{};
		NSpMat newtonHessianElasticity{};
		NSpMat newtonHessianBending{};
		NSpMat newtonHessianCollision{};
		NSpMat newtonHessianAll{};
		NVecDynamic newtonForce{};
		std::vector<NFloatingType> elasticEnergy{};
		NFloatingType uselessHolder{};
		DirectSolverType solverDirect{};
		CGSolverType solverCG{};

		NVecDynamic Ndx;

		void analyzePattern(bool makeCompressed = false);
		virtual void solve(bool patternChanged, bool handleCollision);

		// configurations
		int solverType{ 0 };

		int cgMaxIterations{ 300 };
		FloatingType cgTolerance{ 1e-7 };
		std::vector<IdType> meshOffsets{};
	};

	struct TetMeshNewtonAssembler : public BaseNewtonAssembler {
		void initialize(std::vector<TetMeshFEM::SharedPtr> meshes);

		std::vector<NMat12> elasticHessian{};
		std::vector<NVec12> elasticForce{};

		size_t numAllTets{};
	};

	struct TriMeshCollisionInfoForNewton {
		Vec3 normal;
		Vec4 barys;
		FloatingType energy, lambda, d2E_dDdD;
		bool contacting;
	};

	struct TriMeshNewtonAssembler : public BaseNewtonAssembler {
		// only need to be called once per simulation, it initializes the elastic Hessian and force
		// solverType: 0 for direct solverDirect, 1 for CG solverDirect
		void initialize(std::vector<TriMeshFEM::SharedPtr> meshes_in, int solverType_in = 0);
		// need to be called after each collision detection, it updates the collision Hessian and force
		//void analyzeCollision(std::vector<TriMeshFEM::SharedPtr> meshes);

		void analyzeCollision(const std::vector<std::vector<ClothVFContactQueryResult>>& vfCollisions,
			const std::vector<std::vector<ClothEEContactQueryResult>>& eeCollisions);

		void updatePositions(const VecDynamic& dx);

		void computeBendingHessian();
		std::vector<NMat9> elasticHessian{};
		std::vector<NVec9> elasticForce{};
		std::vector<NVec12> bendingForce{};

		std::vector<NTriplet> newtonHessianTripletsVFCollision;
		std::vector<NTriplet> newtonHessianTripletsEECollision;

		std::vector<std::vector<std::vector<TriMeshCollisionInfoForNewton>>> vfCollisionInfos{};
		std::vector<std::vector<std::vector<TriMeshCollisionInfoForNewton>>> eeCollisionInfos{};

		// pointers to the diagonal blocks of the Hessian matrix
		// nMesh x nVertices x (nContact * 12 * 12)
		// std::vector<std::vector<std::vector<NFloatingType*>>> vfCollisionHessianBlockPtrs{};

		std::vector<TriMeshFEM::SharedPtr> meshes{};

		size_t numAllTris{};


	};


	template<typename T>
	inline bool isNull(const Eigen::SparseMatrix<T>& mat, int row, int col)
	{
		for (Eigen::SparseMatrix<T>::InnerIterator it(mat, col); it; ++it) {
			if (it.row() == row) return false;
		}
		return true;
	}


	inline void addedToSparse3x3Block(NSpMat& m, IdType rowStart, IdType colStart, const Mat3 h) {
		assert(rowStart >= 0 && colStart >= 0 && rowStart + 3 <= m.rows() && colStart + 3 <= m.cols());
		for (IdType i = 0; i < 3; i++) {
			for (IdType j = 0; j < 3; j++) {
				assert(!isNull(m, rowStart + i, colStart + j));
				//RELEASE_ASSERT(!isNull(m, rowStart + i, colStart + j));
				m.coeffRef(rowStart + i, colStart + j) += h(i, j);
			}
		}
	}

}
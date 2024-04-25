#pragma once

#include "../Types/Types.h"
#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"

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

#ifdef USE_MKL
	using DirectSolverType = Eigen::PardisoLDLT<NSpMat>;
	using CGSolverType = Eigen::ConjugateGradient<NSpMat, Eigen::Lower | Eigen::Upper>;
#else
	using DirectSolverType = Eigen::SimplicialLDLT<NSpMat>;
	using CGSolverType = Eigen::ConjugateGradient<NSpMat, Eigen::Lower | Eigen::Upper>;
#endif

	struct BaseNewtonAssembler {

		size_t numAllVertices{};
		size_t numAllEdges{};	

		std::vector<NTriplet> newtonHessianTripletsElasticity;

		// pointers to the diagonal blocks of the Hessian matrix
		std::vector<std::vector<NFloatingType*>> diagonalHessianBlockPtrs{};
		// pointers to the off-diagonal blocks of the Hessian matrix, not necessarily corresponding to the edges
		std::vector<std::vector<NFloatingType*>> offDiagonalHessianBlockPtrs{};

		NFloatingType newtonEnergy{};
		NSpMat newtonHessian{};
		NVecDynamic newtonForce{};
		std::vector<NFloatingType> elasticEnergy{};
		NFloatingType uselessHolder{};
		DirectSolverType solverDirect{};
		CGSolverType solverCG{};

		NVecDynamic Ndx;

		virtual void solve();

		// configurations
		int solverType{ 0 };

		int cgMaxIterations{ 300 };
		FloatingType cgTolerance{ 1e-7 };
	};

	struct TetMeshNewtonAssembler : public BaseNewtonAssembler {
		void initialize(std::vector<TetMeshFEM::SharedPtr> meshes);
		
		std::vector<NMat12> elasticHessian{};
		std::vector<NVec12> elasticForce{};

		size_t numAllTets{};
	};

	struct TriMeshNewtonAssembler : public BaseNewtonAssembler {
		// only need to be called once per simulation, it initializes the elastic Hessian and force
		// solverType: 0 for direct solverDirect, 1 for CG solverDirect
		void initialize(std::vector<TriMeshFEM::SharedPtr> meshes_in, int solverType_in=0);
		// need to be called after each collision detection, it updates the collision Hessian and force
		//void analyzeCollision(std::vector<TriMeshFEM::SharedPtr> meshes);

		// after calling initialize() and analyzeCollision(), call this function to make the Hessian matrix
		void makeHessian(bool makeCompressed=false);

		void updatePositions(const VecDynamic& dx);

		std::vector<NMat9> elasticHessian{};
		std::vector<NVec9> elasticForce{};
		std::vector<NVec12> bendingForce{};

		std::vector<NTriplet> newtonHessianTripletsCollision;
		std::vector<std::vector<NFloatingType*>> collisionHessianBlockPtrs{};

		std::vector<TriMeshFEM::SharedPtr> meshes{};

		size_t numAllTris{};


	};

	

}
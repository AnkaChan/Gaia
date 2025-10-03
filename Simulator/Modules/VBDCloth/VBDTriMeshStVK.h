#pragma once

#include <TriMesh/TriMesh.h>
#include <Framework/BaseClothSimPhsicsFramework.h>
#include "VBDClothPhysicsParameters.h"

#include "./SpatialQuery/MeshClosestPointQuery.h"

#include "VBDBaseTriMesh.h"

namespace GAIA {

	struct VBDObjectParamsTriMeshStVk : public VBDObjectParamsTriMesh {
		typedef std::shared_ptr<VBDObjectParamsTriMeshStVk> SharedPtr;
		typedef VBDObjectParamsTriMeshStVk* Ptr;

		FloatingType miu = 1e4f;
		FloatingType lambda = 1e4f;
		FloatingType dampingStVK = 0.f;
		FloatingType dampingBending = 0.f;

		// 0: inertia, 1: previous position; 
		// 1: adaptive
		// 2: inertia without external force, 
		// 3: symplectic Euler
		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};

	struct VBDTriMeshStVk : public VBDBaseTriMesh
	{
		typedef std::shared_ptr<VBDTriMeshStVk> SharedPtr;
		typedef VBDTriMeshStVk* Ptr;

		void initialize(TriMeshParams::SharedPtr inObjectParams, BaseClothPhsicsFramework::Ptr in_pPhysics);

		//void solve();
		void solverIteration();

		Vec3 greenStrainVoigt(const Mat3x2& F);
		Mat2 greenStrain(const Mat3x2& F);
		void intializeEdgeQuadraticForms();

		// virtual void evaluateInternalForce();
		// void evaluateInternalForcePerTriangle();
		// void evaluateInternalForcePerVertex();
		// void evaluateInternalForcePerEdge_bendingAngle();

		virtual void computeElasticityForceHessianForFace(int tetId, FloatingType& energy, Vec9& force, Mat9& hessian, bool psd = false) override;
		virtual void computeElasticityEnergyForFace(int iFace, FloatingType& energy) override;

		// quardratic bending
		// virtual void evaluateInternalForcePerEdge_quadratic();

		virtual void VBDStep(int iV);

		virtual FloatingType evaluateVertexMeritEnergy(int iV, FloatingType& meInertia, FloatingType& meElastic_stvk, FloatingType& meElastic_bending);
		virtual void accumlateMaterialForceAndHessian(int iV, Vec3& force, Mat3& hessian);

		void accumlateStVKForceAndHessian(int iV, Vec3& force, Mat3& hessian);
		void accumlateBendingForceAndHessian(int iV, Vec3& force, Mat3& hessian);

		FloatingType backTracingLineSearchVBD(IdType vId, const Vec3& dx, FloatingType E0, FloatingType alpha, FloatingType c, FloatingType tau,
			int maxNumIters);

		FloatingType evaluatedAvgGradNorm();
		FloatingType evaluatedRegistrationEnergy();
		FloatingType evaluatedMeritEnergyElastic(FloatingType& meElastic_stvk, FloatingType& meElastic_bending);
		const VBDObjectParamsTriMeshStVk& objectParams() { return *(VBDObjectParamsTriMeshStVk*)(pObjectParams.get()); };

		void clearGradient();
		// void accumulateMaterialGradient(const bool withInertia);

		int tearMesh(IdType v1, IdType v2) override;

		// size: 4 * numFace(),
		// 4*iFace ~ 4*IFace + 3 entries store the diagonal blocks of the Hessian matrix of the energy wrt. the edges's 4 vertices
		// not needed for quadratic bending
		std::vector<Mat4> edgeLaplacianQuadraticForms;

		std::vector<TriMeshClosestPointQueryResult> targetPointsQueyResults;

	};
}
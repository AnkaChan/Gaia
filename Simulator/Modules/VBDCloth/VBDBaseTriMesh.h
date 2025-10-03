#pragma once
#include <TriMesh/TriMesh.h>
#include "VBDClothPhysicsParameters.h"
#include "./CollisionDetector/ClothContactDetector.h"
#include "VBDTriMeshConstraints.h"

namespace GAIA {
	struct VBDObjectParamsTriMesh : public TriMeshParams {
		typedef std::shared_ptr<VBDObjectParamsTriMesh> BaseSharedPtr;
		typedef VBDObjectParamsTriMesh* BasePtr;

		FloatingType exponentialVelDamping = 1.f;
		FloatingType constantVelDamping = 0.f;

		// 0: full inertia | 1: adaptive | 2: previous position | 3: inertia within conservative bounds
		int initializationType = 2;

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};

	struct VBDBaseTriMesh : public TriMeshFEM
	{
		typedef std::shared_ptr<VBDBaseTriMesh> BaseSharedPtr;
		typedef VBDBaseTriMesh* BasePtr;

		virtual void initialize(TriMeshParams::SharedPtr inObjectParams, BaseClothPhsicsFramework::Ptr in_pPhysics);

		//virtual void evaluateInternalForce() = 0;
		//virtual void evaluateInternalForcePerEdge_quadratic() = 0;
		virtual void evaluateExternalForce();

		virtual FloatingType evaluateVertexMeritEnergy(int iV, FloatingType& meInertia, FloatingType& meElastic_stvk, FloatingType& meElastic_bending) = 0;
		virtual FloatingType evaluatedMeritEnergy(FloatingType& meInertia, FloatingType& meElastic_stvk, FloatingType& meElastic_bending);
		virtual FloatingType evaluatedMeritEnergyInertia();
		virtual FloatingType evaluatedMeritEnergyElastic(FloatingType& meElastic_stvk, FloatingType& meElastic_bending) = 0;
		// virtual void accumulateMaterialGradient(const bool withInertia) = 0;

		virtual void computeElasticityForceHessianForFace(int tetId, FloatingType& energy, Vec9& force, Mat9& hessian, bool psd = false) = 0;
		virtual void computeElasticityEnergyForFace(int iFace, FloatingType& energy) = 0;

		const VBDObjectParamsTriMesh& objectParamsBase() { return *(VBDObjectParamsTriMesh*)(pObjectParams.get()); };

		const EdgeInfo& getEdgeInfo(const IdType edgeId) const;
		Vec3 edgeLerp(const IdType edgeId, const FloatingType t) const;

		virtual void accumlateMaterialForceAndHessian(int iV, Vec3& grad, Mat3& hessian) = 0;

		virtual void applyInitialStep();

		//void forwardStep_symplecticEuler();
		void forwardStep_intertia();
		void forwardStep_IE_GD();
		//void forwardStep_IE_VBD();

		virtual void VBDStep(int iV) = 0;
		void accumlateInertiaForceAndHessian(int iV, Vec3& grad, Mat3& hessian);
		void accumlateSewingForceAndHessian(int iV, Vec3& force, Mat3& hessian);

		Vec3Block vertexPositionNext(int iV) { return positionsNext.block<3, 1>(0, iV); }
		ConstVec3Block vertexPositionNext(int iV) const { return positionsNext.block<3, 1>(0, iV); }

		virtual void handleInternalForceConstraint();
		virtual void handleVeclocityConstraint();
		void GAIA::VBDBaseTriMesh::clearGradient();

		BaseClothPhsicsFramework::Ptr pPhysics;
		VBDClothPhysicsParameters::Ptr pPhysicsParams;
		
		VecDynamicBool activeCollisionMask;

		TVerticesMat positionsNext;
		TVerticesMat vertexInternalForces;
		TVerticesMat vertexExternalForces;
		TVerticesMat inertia;
		TVerticesMat gradient;
		TVerticesMat acceletration;

		TVerticesMat positionPrevIter;
		TVerticesMat positionPrevPrevIter;

		VecDynamic vertexConvervativeBounds;
		TVerticesMat positionsAtPrevCollisionDetection;

		VecDynamicBool edgeDegenrateMask;

		VecDynamic meritEnergy_PerVertex;
		VecDynamic meritEnergy_Inertia_PerVertex;
		VecDynamic meritEnergy_StVK_PerVertex;
		VecDynamic meritEnergy_Bending_PerVertex;

		VecDynamic meritEnergyElasticPerFace;
		VecDynamic meritEnergyElasticPerEdge;

		// infos related to face-based energy
 		// shape: 9 x numFace()
 		Mat9x triangleInternalForces;
 
 		// infos related to edge-based energy
 		Mat12x edgeInternalForces; // from top to down: eV12Next, eV1, eV2, eV21Next 
 


		// regsitration data


		// contraints
		std::vector<VertexConstraints> vertexConstraints;

		bool hasVelocitiesPrev = false;
		bool hasApproxAcceleration = false;
	};


	inline Vec3 VBDBaseTriMesh::edgeLerp(const IdType edgeId, const FloatingType t) const
	{
		const EdgeInfo& eInfo = getEdgeInfo(edgeId);
		return (1.f - t) * vertex(eInfo.eV1) + t * vertex(eInfo.eV2);
	}

	inline const EdgeInfo& GAIA::VBDBaseTriMesh::getEdgeInfo(const IdType edgeId) const
	{
		return pTopology->edgeInfos[edgeId];
	}

	inline void assembleMembraneForceAndHessian(const Vec6& dE_dF, const Mat6& d2E_dF_dF, CFloatingType m1, CFloatingType m2, Vec3& dE_dx, Mat3& hessian) {

		Eigen::Matrix<FloatingType, 3, 6> HL;
		CFloatingType A1 = dE_dF(0);
		CFloatingType A2 = dE_dF(1);
		CFloatingType A3 = dE_dF(2);
		CFloatingType A4 = dE_dF(3);
		CFloatingType A5 = dE_dF(4);
		CFloatingType A6 = dE_dF(5);

		dE_dx << A1 * m1 + A4 * m2,
			A2* m1 + A5 * m2,
			A3* m1 + A6 * m2;

		HL.row(0) = d2E_dF_dF.row(0) * m1 + d2E_dF_dF.row(3) * m2;
		HL.row(1) = d2E_dF_dF.row(1) * m1 + d2E_dF_dF.row(4) * m2;
		HL.row(2) = d2E_dF_dF.row(2) * m1 + d2E_dF_dF.row(5) * m2;

		hessian.col(0) = HL.col(0) * m1 + HL.col(3) * m2;
		hessian.col(1) = HL.col(1) * m1 + HL.col(4) * m2;
		hessian.col(2) = HL.col(2) * m1 + HL.col(5) * m2;
	}

	inline void assembleMembraneHessian(const Mat6& d2E_dF_dF, CFloatingType m1, CFloatingType m2, Mat3& hessian) {

		Eigen::Matrix<FloatingType, 3, 6> HL;

		HL.row(0) = d2E_dF_dF.row(0) * m1 + d2E_dF_dF.row(3) * m2;
		HL.row(1) = d2E_dF_dF.row(1) * m1 + d2E_dF_dF.row(4) * m2;
		HL.row(2) = d2E_dF_dF.row(2) * m1 + d2E_dF_dF.row(5) * m2;

		hessian.col(0) = HL.col(0) * m1 + HL.col(3) * m2;
		hessian.col(1) = HL.col(1) * m1 + HL.col(4) * m2;
		hessian.col(2) = HL.col(2) * m1 + HL.col(5) * m2;
	}

	inline void assembleMembraneHessianForWholeFace(const Mat6& d2E_dF_dF, const Vec6& ms, Mat9& hessian) {

		Eigen::Matrix<FloatingType, 3, 6> HL;

		for (int i = 0; i < 3; i++)
		{
			for (int j = 0; j < 3; j++)
			{
				CFloatingType m1_i = ms(2 * i);
				CFloatingType m2_i = ms(2 * i + 1);

				CFloatingType m1_j = ms(2 * j);
				CFloatingType m2_j = ms(2 * j + 1);

				auto h_v = hessian.block<3, 3>(3 * i, 3 * j);

				HL.row(0) = d2E_dF_dF.row(0) * m1_i + d2E_dF_dF.row(3) * m2_i;
				HL.row(1) = d2E_dF_dF.row(1) * m1_i + d2E_dF_dF.row(4) * m2_i;
				HL.row(2) = d2E_dF_dF.row(2) * m1_i + d2E_dF_dF.row(5) * m2_i;

				h_v.col(0) = HL.col(0) * m1_j + HL.col(3) * m2_j;
				h_v.col(1) = HL.col(1) * m1_j + HL.col(4) * m2_j;
				h_v.col(2) = HL.col(2) * m1_j + HL.col(5) * m2_j;
			}
		}
	}
} // namespace GAIA
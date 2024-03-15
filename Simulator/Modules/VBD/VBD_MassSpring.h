#pragma once

#include "../Types/Types.h"
#include "VBD_BaseMaterial.h"
#include <memory>

// #include "../ImplicitEuler/HessianAssebler.h"

namespace GAIA{

	//struct ObjectParametersEBDMassSpring : public ObjectParamsVBD
	//{
	//	typedef std::shared_ptr<ObjectParametersEBDMassSpring> SharedPtr;
	//	typedef ObjectParametersEBDMassSpring* Ptr;

	//	ObjectParametersEBDMassSpring() {
	//		materialType = MassSpring;
	//	}

	//	FloatingType springStiffness = 1e3;
	//	FloatingType damping = 0;

	//	bool lengthBasedMass = true;


	//	virtual bool fromJson(nlohmann::json& objectJsonParams);
	//	virtual bool toJson(nlohmann::json& objectJsonParams);

	//};

	//class EBDTetMesh_MassSpring : public VBDBaseTetMesh {
	//public:
	//	typedef std::shared_ptr<EBDTetMesh_MassSpring> SharedPtr;
	//	typedef EBDTetMesh_MassSpring* Ptr;

	//	EBDTetMesh_MassSpring()
	//	{

	//	}

	//	virtual void solve();

	//	virtual void solve_vertexBlockDescent();
	//	virtual void solve_implicitEuler();
	//	virtual void solve_implicitEuler_iteration();

	//	void solve_energyProjection();
	//	void solve_localEnergyConsistency();

	//	void energyProjectionLocal();
	//	void energyProjectionGlobal();
	//	void energyProjectionLocalForVertex(int iV);
	//	void localMeritEnergyVertexBlockCoordinateDesent(int iV);

	//	FloatingType localMeritEnergyVertexBlockCoordinateBackTrackingLineSearch(int iV, FloatingType E0, Vec3 dx, FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters);


	//	virtual void initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF,
	//		VBDPhysicsParameters::SharedPtr inPhysicsParams, BasePhysicFramework* in_pPhysicsFramework);

	//	virtual void rebalanceMassBasedOnSpringLength();

	//	virtual void forwardStepForwardEuler();

	//	virtual FloatingType evaluateInternalPotentialEnergy(Eigen::Ref<VecDynamic> position);
	//	virtual FloatingType evaluateInternalPotentialEnergy(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> perElementEnergy);

	//	virtual FloatingType evaluateMeritEnergyPervertex(int iV);

	//	virtual FloatingType evaluateLocalInternalEnergyPerVertex(int vId);
	//	virtual void evaluateLocalInternalEnergyHessianPerVertex(int vId, Mat3& h);
	//	virtual FloatingType evaluateLocalInternalEnergyPerVertex(int vId, const Vec3& newPos, const Eigen::Ref<VecDynamic> allVerts);
	//	virtual void evaluateInternalForcePerVertex(int vId, Vec3& force, bool addGravity = true);
	//	virtual void evaluateInternalForcePerVertex(int vId, const Vec3& newPos, const Eigen::Ref<VecDynamic> allVerts, Vec3& force, bool addGravity = true);

	//	virtual void evaluateInternalForce(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> internalForceOut,
	//		bool addGravity=true);
	//	void evaluateSpringEnergyPerSpring(Eigen::Ref<VecDynamic> position, Eigen::Ref<VecDynamic> spingEnergy);
	//	
	//	FloatingType evaluateSpringEnergy(Eigen::Ref<VecDynamic> position);

	//	virtual void computeElasiticForce(const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& p4,
	//		const Mat3& DsInv, FloatingType restVol, Vec3& p1Grad, Vec3& p2Grad, Vec3& p3Grad, Vec3& p4Grad);

	//	virtual FloatingType evaluateEnergy(int iTet);

	//	virtual FloatingType evaluateElasticEnergy();

	//	FloatingType adjustedSpringStiffness(int iSpring);

	//	ObjectParametersEBDMassSpring::SharedPtr pObjectParamsMaterial;

	//	HessianAssembler hessianAssembler;


	//public:
	//	std::vector<std::vector<int>> vertexNeighborEdges;
	//	VecDynamic orgLengths;
	//	// VecDynamic springEnergy; // for temporal storage use

	//	FloatingType avgSpringLength;

	//	Mat2xI edges;
	//	size_t nEdges;

	//	const FloatingType massDivisionRatio = 0.5f;
	//};
}
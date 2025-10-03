#pragma once

#include <Framework/BaseClothSimPhsicsFramework.h>
#include <TriMesh/TriMesh.h>
#include "VBDClothPhysicsParameters.h"
#include "VBDTriMeshStVK.h"
#include <SpatialQuery/MeshClosestPointQuery.h>

#include <SolverUtils/NewtonAssembler.h>
#include <SolverUtils/ChebysevAccelerator.h>

namespace GAIA {

	struct ObjectParamsListVBDCloth : public ObjectParamsList {
		virtual ObjectParams::SharedPtr createObjectParam(const std::string& materialName)
		{
			ObjectParams::SharedPtr pObjParams;
			if (materialName == "StVK_triMesh")
			{
				pObjParams = std::make_shared<VBDObjectParamsTriMeshStVk>();

			}
			//else if (materialName == "NeoHookean") {
			//	pObjParams = std::make_shared<ObjectParametersEBDNeoHookean>();
			//}
			else
			{
				std::cout << "Warning!!! Material name: " << materialName << " not recognized! Skipping this material!\n";
				/*std::cout << "Using default material NeoHookean instead!\n";
				objectParams.push_back(std::make_shared<ObjectParamsPBDNeoHookean>());*/
			}

			return pObjParams;
		}
	};

	class VBDClothBaseDeformer;

	struct VBDClothRuntimeStatistics : public RunningTimeStatistics {
		virtual void setToZero() {
			RunningTimeStatistics::setToZero();
			numberCollisionDetection = 0;
			avgForceResiduals.clear();
			maxForceResiduals.clear();
		}
		virtual std::string customString() {
			std::stringstream ss;
			ss << "-----Number of Collision Detection: " << numberCollisionDetection << "\n";
			return ss.str();
		}

		virtual void customToJson(nlohmann::json& j) {
			PUT_TO_JSON(j, numberCollisionDetection);
			PUT_TO_JSON(j, avgForceResiduals);
			PUT_TO_JSON(j, maxForceResiduals);
		}

		int numberCollisionDetection = 0;

		// substeps x iterations
		std::vector<std::vector<FloatingType>> avgForceResiduals;
		std::vector<std::vector<FloatingType>> maxForceResiduals;
	};


	struct VBDClothSimulationFramework : public BaseClothPhsicsFramework
	{
		VBDClothSimulationFramework()
		{
			simulatorName = "VBDClothSim";
		}

		void initialize() override;
		void initializeNewton();
		void initializeParallelGroups();

		IdType findSmallestParallelGroup(std::vector<IdType>& availableColors, const std::vector<std::vector<IdType>>& vertexParallelGroups,
			bool removeSmallestAvailableColor);
		virtual void parseRunningParameters() override;
		void loadConstraints();

		//void intersectionDetection();
		//void initializeIntersectionDetector();


		virtual TriMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, BasePhysicsParams::SharedPtr physicsParaemters, BaseClothPhsicsFramework* pPhysics);

		virtual std::shared_ptr<ObjectParamsList> createObjectParamsList();
		virtual std::shared_ptr<BasePhysicsParams> createPhysicsParams();
		virtual std::shared_ptr<RunningTimeStatistics> createRunningTimeStatistics();
		virtual void runStepSequential();
		virtual void runStep();
		virtual void runStep_CDEveryColor();
		virtual void runStep_CDEveryIter();
		virtual void runStep_Newton();
		// virtual void runStep_blockJacobi();

		void updateBVH(int groupId = 0);
		void applyCollisionDetection();
		void initializeActiveCollisionMask();

		void applyAcceleration(CFloatingType rho);

		void VBDStepWithCollisionDetection(IdType iMesh, IdType vId, bool updateCollisions);
		void VBDStepWithExistingCollisions(IdType iMesh, IdType vId);
		void evaluateVertexForceAndHessian(IdType iMesh, IdType vId, Vec3& force, Mat3& hessian);
		// return val: true if truncation happened 
		bool applyConservativeBoundTruncation(IdType iMesh, IdType vId, Vec3& newPos);

		CFloatingType accumulateVFContactForceAndHessian(IdType contactId, IdType contactVertexOrder, ClothVFContactQueryResult* pCollisionResult,
			bool updateContactInfo, bool apply_friction, Vec3& force, Mat3& hessian);
		CFloatingType accumulateEEContactForceAndHessian(IdType contactId, IdType contactVertexOrder, ClothEEContactQueryResult* pCollisionResult,
			bool updateContactInfo, bool apply_friction, Vec3& force, Mat3& hessian);
		void accumlateBoundaryForceAndHessian(VBDBaseTriMesh* pMesh, IdType iMesh, IdType vId, Vec3& force, Mat3& hessian);

		void accumulateConstraintForceAndHessian(IdType iMesh, IdType vId, Vec3& force, Mat3& hessian);

		void computeFriction(CFloatingType mu, CFloatingType lambda, const Mat3x2& T, const Vec2& u, CFloatingType epsU, Vec3& force, Mat3& hessian);
		void computeContactRepulsiveForce(CFloatingType dis, CFloatingType thickness, FloatingType& dEdD, FloatingType& d2E_dDdD);

		void computeContactEnergy(CFloatingType dis, CFloatingType thickness, FloatingType& E);

		virtual void tearMesh(IdType iMesh, IdType v1, IdType v2);

		virtual void applyInitialGuess();
		virtual FloatingType evaluateEnergy(ConvergenceStats& energyStats);
		virtual FloatingType evaluateCollisionEnergy(ConvergenceStats& energyStats);
		void clearGradient();

		void evaluateConvergence();

		// void preconditionGradients();
		// void accumulateMaterialGradient();
		void accumulateRegistrationGradient();
		void accumulateCollisionGradient();
		void recomputeContacts();
		void queryClosestTargetPoint();
		// void applyGradientDescent(FloatingType step);
		void revertGradientDescent(FloatingType step);
		void printStats(const ConvergenceStats& energyStats, FloatingType energyAll);
		void saveStats(const std::string outFile, std::vector<ConvergenceStats> convergenceStats);
		void updateVelocity();
		bool determineConvergence(const ConvergenceStats& energyStats, FloatingType energyAll_prev,
			FloatingType avgGradNorm_prev, FloatingType currentStepSize);
		virtual bool writeSimulationParameters(nlohmann::json& outPhysicsParams);

		void collisionDetection();
		void applyDeformers();

		size_t numClothes();

		// FloatingType backTracingLineSearchVBD(FloatingType E0,
		// 	FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters, FloatingType minStepSizeGD, ConvergenceStats& stats);

		void saveIntermediateResults();

		void loadTargetMesh(int iFrame);

		// iMesh < baseTriMeshesForSimulation.size()
		VBDBaseTriMesh& getSimulatedMesh(int iMesh) {
			return *(VBDBaseTriMesh*)(baseTriMeshesForSimulation[iMesh].get());
		}

		// iMesh < triMeshesAll.size()
		// this includes the collider mesh
		TriMeshFEM& getMesh(int iMesh) {
			return *(TriMeshFEM*)(triMeshesAll[iMesh].get());
		}

		VBDObjectParamsTriMesh& getObjectParam(int iObj)
		{
			return objectParamsList->getObjectParamAs<VBDObjectParamsTriMesh>(iObj);
		}
		VBDClothPhysicsParameters& physicsParams()
		{
			return *(VBDClothPhysicsParameters*)(basePhysicsParams.get());
		}
		VBDClothRuntimeStatistics& runtimeStatistics()
		{
			return *(VBDClothRuntimeStatistics*)baseTimeStatistics.get();
		}
		// data
	public:

		std::vector<TriMeshFEM::SharedPtr> targetMeshes;
		MeshClosestPointQueryParameters::SharedPtr pTargetMatcherParameter;
		MeshClosestPointQuery::SharedPtr pTargetMatcher;

		std::vector<std::shared_ptr<VBDClothBaseDeformer>> deformers;

		// for temporally storing the original positions in the process of back tracking line search
		std::vector<TVerticesMat> orgPosistions;
		// convergenceStats for a frame
		// convergenceStats.size() = iterations
		std::vector<ConvergenceStats> convergenceStats;

		// collision data

		ClothVFContactQueryResult& getVFContactResult(IdType iMesh, IdType iV)
		{
			return clothClothVFContactQueyResults[iMesh][iV];
		}
		ClothEEContactQueryResult& getEEContactResult(IdType iMesh, IdType iE)
		{
			return clothClothEEContactQueyResults[iMesh][iE];
		}
		ClothVFContactQueryResult& getFVContactResult(IdType iMesh, IdType faceId)
		{
			return clothClothFVContactQueyResults[iMesh][faceId];
		}
		std::vector<std::vector<ClothVFContactQueryResult>> clothClothVFContactQueyResults;
		std::vector<std::vector<ClothVFContactQueryResult>> clothClothFVContactQueyResults;
		std::vector<std::vector<ClothEEContactQueryResult>> clothClothEEContactQueyResults;

		bool contactInfoUpdated = false;
		bool collisionDetectionRequired = true;

		// nGroups x (2 * nVertices)
		// each groups has this structure: iMesh1, iVertex1, iMesh2, iVertex2, ...
		std::vector<std::vector<IdType>> vertexParallelGroups;

	public:
		std::shared_ptr<TriMeshNewtonAssembler> pNewtonAssembler;
		VecDynamic newtonDx;


		void newtonEvaluateElasticForceAndHessian();
		void newtonEvaluateCollisionForceAndHessian();
		NFloatingType newtonEvaluateCollisionEnergy();

		void fillNewtonSystem();
		void fillNewtonForce();
		void fillNewtonHessianDiagonal();
		void fillNewtonHessianOffDiagonal();
		void fillNewtonCollisionForceAndHessian();
		// use pre-computed collision info
		void fillNewtonCollisionForceAndHessianV2();
		void newtonEvaluateElasticEnergy();

		void newtonConservativeStepCulling();

		// return value: true if actively contact
		bool accumulateVFContactForceAndHessianAllVerts(IdType contactId, ClothVFContactQueryResult* pCollisionResult,
			bool updateContactInfo, bool apply_friction, FloatingType& lambda, FloatingType& d2E_dDdD, Vec3& normal, Vec4& barys);
		bool accumulateEEContactForceAndHessianAllVerts(IdType contactId, ClothEEContactQueryResult* pCollisionResult,
			bool updateContactInfo, bool apply_friction, FloatingType& lambda, FloatingType& d2E_dDdD, Vec3& normal, Vec4& barys);

		bool accumulateVFContactEnergy(IdType contactId, ClothVFContactQueryResult* pCollisionResult, FloatingType& energy);
		bool accumulateEEContactEnergy(IdType contactId, ClothEEContactQueryResult* pCollisionResult, FloatingType& energy);


		NFloatingType newtonEvaluateMeritEnergy(NFloatingType& eInertia, NFloatingType& eElastic, NFloatingType& eBending, NFloatingType& eContact, bool elasticReady = false);
		NFloatingType newtonLineSearch(const VecDynamic& dx, NFloatingType E0, FloatingType alpha,
			FloatingType c, FloatingType tau, int maxNumIters, FloatingType& stepSizeOut);
	};
}
#pragma once
#include <MeshFrame/Utility/Parser.h>

#include "../Framework/BasePhysicsFramework.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

#include "VBD_BaseMaterial.h"
#include "VBD_NeoHookean.h"
#include "VBD_MassSpring.h"
#include "VBDPhysicsParameters.h"

#include "ActiveCollisionList.h"
#include "VBDPhysicsTest.h"

#include "VBD_CollisionInfo.h"

#include "../SolverUtils/GD_SolverUtilities.h"
#include "../SolverUtils/LineSearchUtilities.h"
#include "../SolverUtils/NewtonAssembler.h"	
#ifdef USE_MKL
#include <Eigen/PardisoSupport>
#endif


//#define TURN_ON_SINGLE_THREAD_DEBUG
#include "../Parallelization/CPUParallelization.h"


namespace GAIA {
	struct ObjectParamsListVBD : ObjectParamsList {
		virtual ObjectParams::SharedPtr createObjectParam(const std::string& materialName);
	};

	class VBDBaseDeformer;

	struct VBDPhysics : public BasePhysicFramework {
	public:
		typedef std::shared_ptr<VBDBaseTetMesh> BaseMaterialPtr;
		typedef std::shared_ptr<VBDTetMeshNeoHookean> NeoHookeanMaterialPtr;
		//typedef std::shared_ptr<EBDTetMesh_MassSpring> MassSpringMaterialPtr;

		inline ObjectParamsVBD& getObjectParam(int iObj);
		inline RunningTimeStatistics& timeStatistics();
		template<typename MeshType>
		inline std::shared_ptr<MeshType> getTetMeshSharedPtrAs(int iMesh);

		template<typename MeshType>
		inline MeshType& getTetMeshAs(int iMesh);

		void loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder);

		virtual void initialize();
		virtual void initializeNewton();
		virtual void initializeGPU();
		virtual void initializeGPU_cpuDebugData();
		virtual void enableModels();
		virtual void initializeGD();

		void sortVertexColorGroupsByMortonCode();


		virtual TetMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, std::shared_ptr<TetMeshMF> pTMeshMF,
			BasePhysicsParams::SharedPtr physicsParaemters);

		void timeStepEqualize();
		void runStep();
		void runStepGPU();
		void runStepGPU_allInOneSweep();
		void runStepGPU_acceleratedGS();
		void runStepGPU_GD();
		void runStepGPU_debugOnCPU();
		void runStepGPUNoCollision();
		void runStepNewton();

		// record VBD graph
		void recordVBDSolveGraph();

		void applyDeformers();

		void runStep_serialCollisionHandling();
		void runStep_hybridCollisionHandling();

		void prepareCollisionDataCPU();
		void prepareCollisionDataGPU();
		void dcd();
		void intermediateDCD();
		void ccd();
		void intermediateCCD();
		void intermediateCollisionDetection();
		void updateDCDBVH(bool rebuildTetMeshScene, bool rebuildSurfaceScene);
		void updateCCDBVH(bool rebuildScene);
		void updateAllCollisionInfos();
		void updateCollisionInfo(VBDCollisionDetectionResult& collisionResult);
		void solveCollisionsSequentially();

		void VBDStep(TetMeshFEM* pMesh, IdType vertexId);
		void VBDStepWithCollision(TetMeshFEM* pMesh, IdType meshId, IdType vertexId, bool apply_friction = false);

		void updateVelocities();
		void updateVelocitiesGPU();
		void updateVelocity(TetMeshFEM* pMesh, IdType vertexId);
		void dampVelocities();
		// IPC friction
		void computeVertexFriction(CFloatingType mu, CFloatingType lambda, const Mat3x2& T, const Vec2& u, CFloatingType epsU, Vec3& force, Mat3& hessian);

		Vec3 applyFrictinalVelocityDamping(Vec3 velocity, Vec3& contactNormal, CFloatingType vertInvMass,
			CFloatingType frictionRatio, CFloatingType contactForce, CFloatingType dt);

		void accumlateDampingForceAndHessian(TetMeshFEM* pMesh, IdType vertexId, Vec3& force, Mat3& hessian, const Mat3& K);

		void solveBoxBoundaryConstraintForVertex(TetMeshFEM* pMesh, IdType vertexId);
		void accumlateBoundaryForceAndHessian(TetMeshFEM* pMesh, IdType meshId, IdType vertexId, Vec3& force, Mat3& hessian, bool apply_friction = false);
		void applyBoudnaryFriction(TetMeshFEM* pMesh, IdType vertexId);
		void applyCollisionFriction(VBDBaseTetMesh* pMesh, IdType meshId, IdType vertexId);
		// FloatingType computeCollisionRepulsiveForce(VBDBaseTetMesh* pMesh, IdType meshId, IdType vertexId);
		FloatingType computeCollisionRepulsiveForcePerCollision(int collisionId, int collisionVertexOrder,
			VBDCollisionDetectionResult& collisionResult, Vec3& contactNormal);

		// collision Id: 0~3; 0~2 the 3 vertices of the colliding face
		void accumlateCollisionForceAndHessian(TetMeshFEM* pMesh, IdType meshId, IdType vertexId,
			Vec3& force, Mat3& hessian, bool apply_friction = false);
		// update all the related collision infos for vertex
		void updateCollisionInfoForVertex(TetMeshFEM* pMesh, IdType meshId, IdType vertexId);
		void accumlateCollisionForceAndHessianPerCollision(int collisionId, int collisionVertexOrder,
			VBDCollisionDetectionResult& collisionResult, Vec3& force, Mat3& hessian, IdType meshId, IdType vertexId, bool apply_friction = false);

		VBDCollisionDetectionResult& getCollisionDetectionResultFromSurfaceMeshId(int meshId, int vertexSurfaceMeshId);
		VBDCollisionDetectionResult& getCollisionDetectionResultFromTetMeshId(int meshId, int vertexTetMeshId);
		// how a vertex is connected to other points through collision
		CollisionRelationList& getCollisionRelationList(int meshId, int vertexId);

		virtual std::shared_ptr<ObjectParamsList> createObjectParamsList();
		virtual std::shared_ptr<BasePhysicsParams> createPhysicsParams();
		virtual std::shared_ptr<CollisionDetectionParamters> createCollisionParams();
		virtual std::shared_ptr<RunningTimeStatistics> createRunningTimeStatistics();

		VBDPhysicsParameters& physicsParams();
		CollisionDetectionParamters& collisionParams();

		size_t numVertexParallelGroups() { return vertexParallelGroups.size(); }
		size_t sizeVertexParallelGroup(int iGroup) { return vertexParallelGroups[iGroup].size() / 2; }
		int32_t* getVertexParallelGroupsHeadGPU(int iGroup) { return vertexParallelGroupHeadsGPU[iGroup]; }

		size_t numTetParallelGroups() { return tetParallelGroups.size(); }
		size_t sizeTetParallelGroup(int iGroup) { return tetParallelGroups[iGroup].size() / 4; }
		int32_t* getTetParallelGroupsHeadGPU(int iGroup) { return tetParallelGroupHeadsGPU[iGroup]; }

		VBDPhysicsDataGPU* getVBDPhysicsDataGPU() { return pPhysicsDataGPUBuffer->getData(); };

		void syncAllToGPU(bool sync = false);
		void syncAllToCPU(bool sync = false);
		void syncAllToGPUVertPosOnly(bool sync = false);
		void syncAllToCPUVertPosOnly(bool sync = false);

		// test code
		void testGPUCollisionHandlingCode();
		void prepareCollisionDataCPUAndGPU();
		void printCPUCollisionDataForVertex(IdType meshId, IdType vertexId);
		void printCPUCollisionData(VBDCollisionDetectionResult& colResult, int iIntersection = 0);
		void compareCPUandGPU(Vec3& forceCPU, Mat3& hessianCPU, Vec3& forceGPU, Mat3& hessianGPU);

		// debug
		void outputPosVel();
		void outputForces();
		void clearForces();
		void evaluateConvergence();
		void evaluateConvergenceGPU();

		// accelerator
		FloatingType getAcceleratorOmega(int order, CFloatingType pho, CFloatingType prevOmega);
		void recordInitialPositionForAccelerator(bool sync);

		// data
	public:
		// tetmesh and corresponding materials
		std::vector<VBDBaseTetMesh::SharedPtr> tMeshes;
		std::vector<ObjectParamsVBD::SharedPtr> objectParamsVBD;
		VBDPhysicsParameters::SharedPtr physicsParamsVBD;

		// debug
		std::vector<TVerticesMat> MaterialForce{};
		std::vector<TVerticesMat> InertiaForce{};
		std::vector<TVerticesMat> boundaryFrictionForce{};
		std::vector<TVerticesMat> boundaryCollisionForce{};
		std::vector<TVerticesMat> frictionForce{};
		std::vector<TVerticesMat> collisionForce{};

		std::vector<std::shared_ptr<VBDBaseDeformer>> deformers;

		// GPU data
		std::vector<VBDBaseTetMeshGPU*> tMeshesGPU;
		ManagedBuffer<VBDBaseTetMeshGPU*>::SharedPtr tetMehesGPUBuffer;
		VBDPhysicsDataGPU vbdPhysicsDataGPUCPUBuffer;
		std::shared_ptr<DeviceClassBuffer<VBDPhysicsDataGPU>> pPhysicsDataGPUBuffer;

		std::vector<ManagedBuffer<int32_t>::SharedPtr> vertexParallelGroupsBuffer;
		std::vector<int32_t*> vertexParallelGroupHeadsGPU;

		std::vector<ManagedBuffer<int32_t>::SharedPtr> tetParallelGroupsBuffer;
		std::vector<int32_t*> tetParallelGroupHeadsGPU;

		// for the operation that allows full parallelization by vertex, like updating the velocity
		// meshId1, vertexId1, meshId2, vertexId2, meshId3, vertexId3 ...
		ManagedBuffer<int32_t>::SharedPtr vertexAllParallelGroupsBuffer;

		std::vector<ManagedBuffer<int32_t>::SharedPtr> activeCollisionsEachParallelGroupBuffer;

		// for debugging GPU
		std::vector<VBDBaseTetMeshGPU*> tetMehesGPUBuffer_forCPUDebug;
		VBDPhysicsDataGPU vbdPhysicsDataGPU_forCPUDebug;

		// for GPU executions
		cudaStream_t cudaStream;
		bool graphCreated = false;
		cudaGraph_t VBDSolveGraph;
		cudaGraphExec_t VBDSolveInstance;
		// nMesh x nVertices
		std::vector<std::vector<VBDCollisionDetectionResult>> collisionResultsAll;

		// nGroups x (2 * nVertices)
		// each groups has this structure: iMesh1, iVertex1, iMesh2, iVertex2, ...
		std::vector<std::vector<IdType>> vertexParallelGroups;
		// nGroups x (4 * nTets)
		// each groups has this structure: iMesh1, tetId, vertexOrder, vertexId, ...
		std::vector<std::vector<IdType>> tetParallelGroups;

		ActiveCollisionList activeColllisionList;

	private:
		// for line search
		LineSearchUtilities::SharedPtr pLineSearchUtilities;
		FloatingType evaluateMeritEnergyGPU(FloatingType& eInertia, FloatingType& eElastic, bool elasticReady = false);

		TetMeshNewtonAssembler::SharedPtr pNewtonAssembler;

		void computeElasticForceHessian();
		void computeElasticEnergy();
		void fillNewtonSystem();
		void fillNewtonForce();
		void fillNewtonHessianDiagonal();
		void fillNewtonHessianOffDiagonal();
		void updatePositions(const VecDynamic& dx);
		NFloatingType evaluateMeritEnergy(NFloatingType& eInertia, NFloatingType& eElastic, bool elasticReady = false);
		NFloatingType newtonLineSearch(const VecDynamic& dx, NFloatingType E0, FloatingType alpha,
			FloatingType c, FloatingType tau, int maxNumIters, FloatingType& stepSizeOut);

		//For GD solver
		GD_SolverUtilities::SharedPtr pGDSolverUtilities;
		void GDBackupPositions(bool sync, CFloatingType omega);
		void GDRevertToBackupPositions(bool sync, FloatingType& omega);
		FloatingType GDLineSearch(GD_SolverUtilities::SharedPtr pGDSolverUtilities, FloatingType& E_prev,
			FloatingType& EInertia_prev, FloatingType& EElastic_prev,
			FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters);

	};

	inline VBDPhysicsParameters& VBDPhysics::physicsParams()
	{
		return *physicsParamsVBD;
	}

	inline CollisionDetectionParamters& VBDPhysics::collisionParams()
	{
		return *baseCollisionParams;
	}

	inline void VBDPhysics::syncAllToGPU(bool sync)
	{
		for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			tMeshes[iMesh]->syncToGPU(false, cudaStream);
		}
		if (sync)
		{
			CHECK_CUDA_ERROR(cudaDeviceSynchronize());
		}

		// sync collision active list to GPU
		for (int iGroup = 0; iGroup < activeCollisionsEachParallelGroupBuffer.size(); iGroup++)
		{
			size_t numActiveCollisions = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
			activeCollisionsEachParallelGroupBuffer[iGroup]->toGPU(numActiveCollisions * 2, false, cudaStream);
		}

		// test GPU collision data
		// sync collision active list to GPU
		//std::cout << "-----------------------------------------\n" 
		//	<< "Test GPU collision data" << "\n";
		//for (int iGroup = 0; iGroup < activeCollisionsEachParallelGroupBuffer.size(); iGroup++)
		//{
		//	std::cout << "Printing collision data for parallel group:" << iGroup << "\n";
		//	size_t numActiveCollisions = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
		//	printGPUCollisionData(numActiveCollisions, 
		//		activeCollisionsEachParallelGroupBuffer[iGroup]->getGPUBuffer(),
		//		pPhysicsDataGPUBuffer->getData()
		//	);
		//}


	}

	inline void VBDPhysics::syncAllToCPU(bool sync)
	{
		for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			tMeshes[iMesh]->syncToCPU(false, cudaStream);
		}
		if (sync)
		{
			CHECK_CUDA_ERROR(cudaDeviceSynchronize());
		}
	}

	inline void VBDPhysics::syncAllToGPUVertPosOnly(bool sync)
	{
		for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			tMeshes[iMesh]->syncToGPUVertPosOnly(false, cudaStream);
		}
		if (sync)
		{
			CHECK_CUDA_ERROR(cudaDeviceSynchronize());
		}
	}

	inline void VBDPhysics::syncAllToCPUVertPosOnly(bool sync)
	{
		for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			tMeshes[iMesh]->syncToCPUVertPosOnly(false, cudaStream);
		}
		if (sync)
		{
			CHECK_CUDA_ERROR(cudaDeviceSynchronize());
		}
	}



	inline ObjectParamsVBD& VBDPhysics::getObjectParam(int iObj)
	{
		return objectParamsList->getObjectParamAs<ObjectParamsVBD>(iObj);
	}

	inline RunningTimeStatistics& VBDPhysics::timeStatistics()
	{
		// TODO: insert return statement here
		return *(RunningTimeStatistics*)baseTimeStatistics.get();
	}

	template<typename MeshType>
	inline std::shared_ptr<MeshType> VBDPhysics::getTetMeshSharedPtrAs(int iMesh)
	{
		return std::static_pointer_cast<MeshType>(basetetMeshes[iMesh]);
	}

	template<typename MeshType>
	inline MeshType& VBDPhysics::getTetMeshAs(int iMesh)
	{
		*(MeshType*)basetetMeshes[iMesh].get();
	}

}
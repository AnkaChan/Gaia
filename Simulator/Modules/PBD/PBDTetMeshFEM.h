#pragma once
#include "cuda_runtime.h"

#include "../TetMesh/TetMeshFEM.h"
#include "PBDTetMeshFEMGPU.h"
#include "../PBD/PBDPhysicsParameters.h"

namespace GAIA {
	struct PBDPhysics;

	struct ObjectParamsPBD : public ObjectParams {

		typedef std::shared_ptr<ObjectParamsPBD> SharedPtr;
		typedef ObjectParamsPBD* Ptr;

		FloatingType exponentialVelDamping = 1.0;
		FloatingType constantVelDamping = 0.0;

		FloatingType friction_dynamic = 0.12;
		FloatingType friction_static = 0.15;


		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};

	struct PBDTetMeshFEM : public TetMeshFEM
	{
		typedef std::shared_ptr<PBDTetMeshFEM> SharedPtr;
		typedef PBDTetMeshFEM* Ptr;

		virtual void initialize(ObjectParams::SharedPtr inMaterialParams, TetMeshMF::SharedPtr pTM_MF, PBDPhysics* inPPBDPhysics);
		
		PBDPhysics* pPBDPhysics = nullptr;

		virtual void applyInitialGuess();
		virtual void solveMaterialConstraint();
		virtual void solveBoundaryConstraint();
		virtual void solveInversionConstraint();

		virtual void solveMaterialConstraintGPU_ParallelizationGroup(int iGroup, int numThreads, cudaStream_t stream) = 0;
		virtual void solveMaterialGPUAggregated(int numThreads, cudaStream_t stream) = 0;
		virtual void solveMaterialConstraintGPU_debugOnCPU(PBDTetMeshFEMGPU* pExternGPUTMeshOnCPU) = 0;
		virtual size_t numAllParallelizationGroups() = 0;

		virtual size_t numTetParallelizationGroups() = 0;
		virtual size_t getTetParallelizationGroupSize(int32_t iGroup) = 0;
		virtual int32_t* getTetParallelizationGroupBufferGPU(int32_t iGroup) = 0;
		
		virtual void syncToGPU(bool sync = true, cudaStream_t stream = 0) = 0;
		virtual void syncToCPU(bool sync = true, cudaStream_t stream = 0) = 0;
		virtual void syncToCPUVertPosOnly(bool sync = true, cudaStream_t stream = 0) = 0;
		virtual void syncToGPUVertPosOnly(bool sync = true, cudaStream_t stream = 0) = 0;
		virtual void syncToCPUInversionSignOnly(bool sync = true, cudaStream_t stream = 0) = 0;

		virtual FloatingType evaluateMeritEnergy();
		virtual FloatingType evaluateInertiaEnergy();
		virtual FloatingType evaluateElasticEnergy() = 0;

		FloatingType dt;
		TVerticesMat inertia;
	protected:
		// hide this intializer, because PBDTetMeshFEM must be initialized with the pointer to PBDPhysics
		using TetMeshFEM::initialize;
	};

}
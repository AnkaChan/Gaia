#pragma once

#include "PBDTetMeshFEM.h"
#include "PBDTetMeshFEMShared.h"
#include "PBDTetMeshMassSpringCompute.h"

namespace GAIA {
	struct ObjectParamsPBDMassSpring : public ObjectParamsPBD {
		ObjectParamsPBDMassSpring() {
			materialType = MassSpring;
		}

		FloatingType springCompliance = 1e-4f;
		FloatingType springDamping = 0.01f;

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);

		typedef std::shared_ptr<ObjectParamsPBDMassSpring> SharedPtr;
	};

	class PBDTetMeshMassSpring;

	struct PBDTetMeshFEMShared_MassSpring : public PBDTetMeshFEMShared
	{
		std::vector<ManagedBuffer<int32_t>::SharedPtr> edgesColoringEachCategoryBuffer;

		ManagedBuffer<FloatingTypeGPU>::SharedPtr springLambdasBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr orgLengthsBuffer;
		ManagedBuffer<int32_t>::SharedPtr edgesBuffer;

		ManagedBuffer<int32_t*>::SharedPtr edgesColoringCategoriesPointersBuffer;
		ManagedBuffer<int32_t>::SharedPtr edgesColoringEachCategorySizeBuffer;

		virtual void copyAllToGPU();
		virtual void initializeGPUTetMesh(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);
		virtual void initializeGPUTetMeshOnCPU(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);
		// virtual void syncToGPU();
	};

	struct PBDTetMeshMassSpring : public PBDTetMeshFEM, public PBDTetMeshFEMShared_MassSpring
	{
		typedef std::shared_ptr<PBDTetMeshMassSpring> SharedPtr;
		typedef PBDTetMeshMassSpring* Ptr;
		// for the convenience of use
		ObjectParamsPBDMassSpring::SharedPtr pObjectParamsMaterial;

		std::shared_ptr<DeviceClassBuffer<PBDTetMeshFEMGPUMassSpring>> pTetMeshGPUBuffer;
		PBDTetMeshFEMGPUMassSpring* getTetMeshGPU();

		// GPUTetMesh prepared for GPU computation, where each buffer points to the GPU memory
		PBDTetMeshFEMGPUMassSpring tetMeshGPU;
		// GPUTetMesh prepared for CPU computation, where each buffer points to the CPU RAM
		PBDTetMeshFEMGPUMassSpring tetMeshGPU_forCPU;
		//PBDTetMeshFEMGPUMassSpring tetMeshGPU;

		virtual void initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF, PBDPhysics* inPPBDPhysics);

		virtual void applyInitialGuess();

		virtual void syncToGPU(bool sync=true, cudaStream_t stream = 0);
		virtual void syncToCPU(bool sync=true, cudaStream_t stream = 0);
		virtual void syncToCPUVertPosOnly(bool sync=true, cudaStream_t stream = 0);
		virtual void syncToGPUVertPosOnly(bool sync=true, cudaStream_t stream = 0);
		virtual void syncToCPUInversionSignOnly(bool sync = true, cudaStream_t stream = 0);

		// before calling this, call syncToGPU() first
		virtual void solveMaterialGPUAggregated(int numThreads, cudaStream_t stream);
		virtual void solveMaterialConstraintGPU_ParallelizationGroup(int iGroup, int numThreads, cudaStream_t stream);
		virtual void solveMaterialConstraintGPU_debugOnCPU(PBDTetMeshFEMGPU* pExternGPUTMeshOnCPU);

		virtual size_t numAllParallelizationGroups();
		virtual size_t numTetParallelizationGroups();
		virtual size_t getTetParallelizationGroupSize(int32_t iGroup);
		virtual int32_t* getTetParallelizationGroupBufferGPU(int32_t iGroup);

		virtual FloatingType evaluateElasticEnergy();


		void initializeLambda();

		VecDynamic springLambdas;
		VecDynamic orgLengths;

		ObjectParamsPBDMassSpring::SharedPtr getObjectParams();

	};

	inline PBDTetMeshFEMGPUMassSpring* GAIA::PBDTetMeshMassSpring::getTetMeshGPU()
	{
		return pTetMeshGPUBuffer->getData();
	}
}
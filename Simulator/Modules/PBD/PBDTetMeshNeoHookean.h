#pragma once

#include "PBDTetMeshFEM.h"
#include "PBDTetMeshFEMShared.h"
#include "PBDTetMeshNeoHookeanCompute.h"



namespace GAIA {
	struct ObjectParamsPBDNeoHookean : public ObjectParamsPBD {
		ObjectParamsPBDNeoHookean() {
			materialType = NeoHookean;
		}

		FloatingType devCompliance = 1e-4;
		FloatingType volCompliance = 0;

		FloatingType devDamping = 0.01;
		FloatingType volDamping = 0;


		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);

		typedef std::shared_ptr<ObjectParamsPBDNeoHookean> SharedPtr;

	};

	class PBDTetMeshNeoHookean;

	struct PBDTetMeshFEMShared_NeoHookean : public PBDTetMeshFEMShared
	{

		ManagedBuffer<FloatingTypeGPU>::SharedPtr devLambdasBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr volLambdasBuffer;

		FloatingTypeGPU devCompliance;
		FloatingTypeGPU volCompliance;

		virtual void copyAllToGPU();
		virtual void initializeGPUTetMesh(PBDTetMeshFEM* pTetMesh, PBDTetMeshFEMGPU* pTetMeshGPU);
		virtual void initializeGPUTetMeshForDebug(PBDTetMeshFEM* pTetMesh, PBDTetMeshFEMGPU* pTetMeshGPU);

	};

	struct PBDTetMeshNeoHookean : public PBDTetMeshFEM, public PBDTetMeshFEMShared_NeoHookean
	{
		typedef std::shared_ptr<PBDTetMeshNeoHookean> SharedPtr;
		typedef PBDTetMeshNeoHookean* Ptr;
		// buffer of GPU tetMesh
		ObjectParamsPBDNeoHookean::SharedPtr pObjectParamsNeoHookean;
		
		// GPUTetMesh prepared for GPU computation, where each buffer points to the GPU memory
		TetMeshFEMGPU_NeoHookean tetMeshGPU;
		// GPUTetMesh prepared for CPU computation, where each buffer points to the CPU RAM
		// this is for runing GPU code on the cpu
		TetMeshFEMGPU_NeoHookean tetMeshGPU_forCPU;
		
		// inchange of CPU-GPU communication
		std::shared_ptr<DeviceClassBuffer<TetMeshFEMGPU_NeoHookean>> pTetMeshGPUBuffer;
		TetMeshFEMGPU_NeoHookean* getTetMeshGPU();

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
		void initializeLambda();

		virtual FloatingType evaluateElasticEnergy();

		VecDynamic devLambdas;
		VecDynamic volLambdas;

		// for rest stability
		FloatingType lambda;
		FloatingType mu;

		ObjectParamsPBDNeoHookean::SharedPtr getObjectParams();


	};

	inline TetMeshFEMGPU_NeoHookean* GAIA::PBDTetMeshNeoHookean::getTetMeshGPU()
	{
		return pTetMeshGPUBuffer->getData();
	}
}
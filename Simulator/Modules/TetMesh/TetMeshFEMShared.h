#pragma once

#include "TetMeshFEM.h"
#include "TetMeshFEMGPU.h"
#include "CuMatrix/Buffers/ManagedBuffer.h"

#define SET_GPU_MESH_BUFER_PTR(gpuMesh, memberPtr)\
gpuMesh->memberPtr = memberPtr##Buffer->getGPUBuffer();

namespace GAIA {

	struct TetMeshTopologyShared {
		ManagedBuffer<int32_t>::SharedPtr tetVIdsBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborTetsBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborTets_vertexOrderBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborTets_infosBuffer;
		std::shared_ptr<DeviceClassBuffer<TetMeshTopologyGPU>> topologyGPUBuffer;
		std::shared_ptr<TetMeshTopologyGPU> topologyGPUForCPUDebug;

		void initialize(TetMeshTopology* pTopology);
		void initializeForCPUDebug(TetMeshTopology* pTopology);

	};

	// a class manages and syncs data between CPU and GPU
	struct TetMeshFEMShared {

		typedef std::shared_ptr<TetMeshFEMShared> SharedPtr;
		typedef TetMeshFEMShared* Ptr;

		// managed buffers for non-variant data

		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetRestVolumeBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetInvRestVolumeBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertexMassBuffer;

		// flatten vertex data: 3 x nVerts
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPrevPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr velocityBuffer;

		// tet's inversion of restpose transformation: 9 x nTets
		ManagedBuffer<FloatingTypeGPU>::SharedPtr DmInvsBuffer;

		ManagedBuffer<int8_t>::SharedPtr vertexFixedMaskBuffer;

		//// tet's inversion sign buffer:  nTets
		//ManagedBuffer<int8_t>::SharedPtr tetInvertedSignBuffer;
		//// tet's inversion sign buffer:  nVerts
		//ManagedBuffer<int8_t>::SharedPtr verticesInvertedSignBuffer;

		// from the mesh path to tetmesh GPU topologies on GPU memory
		static std::map<std::string, std::shared_ptr<TetMeshTopologyShared>> topologies;
		static std::mutex topologies_lock;
		TetMeshTopologyGPU* getGPUTetMeshTopology(TetMeshFEM* pTetMesh);
		TetMeshTopologyGPU* getCPUTetMeshTopology(TetMeshFEM* pTetMesh);

#ifdef EVAL_TET_ACCESS_COUNT
		ManagedBuffer<int32_t>::SharedPtr tetAccessCountBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertAccessCountBuffer;
#endif //EVAL_TET_ACCESS_COUNT

		// Generate GPU buffers, then copy them to a GPU mesh on CPU, which is ready to be copied to GPU
		void initializeGPUTetMesh(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);
		// Unlike the above one, this generates the GPU TetMesh with CPU buffers
		void initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);


		virtual void copyAllToGPU(bool sync = true, cudaStream_t stream = 0);
		// sync all the mutable arrays, namely vertPos, vertPrevPos, and velocity to GPU
		virtual void syncToGPU(bool sync=true, cudaStream_t stream = 0);
		// sync all the mutable arrays, namely vertPos, vertPrevPos, and velocity to CPU
		virtual void syncToCPU(bool sync = true, cudaStream_t stream = 0);

		void syncToCPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);
		void syncToGPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);

		virtual void syncToCPUInversionSignOnly(bool sync = true, cudaStream_t stream = 0);
	};

	inline void TetMeshFEMShared::initializeGPUTetMesh(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU)
	{
		// std::cout << "Initializing GPU tet mesh for base material.\n";

		TetMeshTopologyGPU* pTopology = getGPUTetMeshTopology(pTetMesh);

		pTetMeshGPU->pTopology = pTopology;

		pTetMeshGPU->activeForCollision = pTetMesh->activeForCollision;
		pTetMeshGPU->activeForMaterialSolve = pTetMesh->activeForMaterialSolve;

		tetRestVolumeBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->tetRestVolume.size(), true, pTetMesh->tetRestVolume.data());
		tetRestVolumeBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetRestVolume);

		tetInvRestVolumeBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->tetInvRestVolume.size(), true, pTetMesh->tetInvRestVolume.data());
		tetInvRestVolumeBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvRestVolume);

		vertexMassBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->vertexMass.size(), true, pTetMesh->vertexMass.data());
		vertexMassBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertexMass);

		DmInvsBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->DmInvs.size(), true, pTetMesh->DmInvs.data());
		DmInvsBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, DmInvs);

		vertPosBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVertPos.size(), true, pTetMesh->mVertPos.data());
		vertPosBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPos);

		vertPrevPosBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVertPrevPos.size(), true, pTetMesh->mVertPrevPos.data());
		vertPrevPosBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPrevPos);

		velocityBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVelocity.size(), true, pTetMesh->mVelocity.data());
		velocityBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, velocity);

		vertexFixedMaskBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->fixedMask.size(),
			true, pTetMesh->fixedMask.data());
		vertexFixedMaskBuffer->toGPU(false);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertexFixedMask);

		/*tetInvertedSignBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->tetsInvertedSign.size(),
			true, pTetMesh->tetsInvertedSign.data());
		tetInvertedSignBuffer->toGPU(false);
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvertedSign);

		verticesInvertedSignBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->verticesInvertedSign.size(),
			true, pTetMesh->verticesInvertedSign.data());
		verticesInvertedSignBuffer->toGPU(false);*/
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, verticesInvertedSign);

#ifdef EVAL_TET_ACCESS_COUNT

		tetAccessCountBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->numTets(), true);
		vertAccessCountBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->numVertices(), true);
		for (size_t iTet = 0; iTet < pTetMesh->numTets(); iTet++)
		{
			tetAccessCountBuffer->getCPUBuffer()[iTet] = 0;
		}
		tetAccessCountBuffer->toGPU();
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetAccessCount);

		for (size_t iVert = 0; iVert < pTetMesh->numVertices(); iVert++)
		{
			vertAccessCountBuffer->getCPUBuffer()[iVert] = 0;
		}
		vertAccessCountBuffer->toGPU();
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertAccessCount);

#endif // !EVAL_TET_ACCESS_COUNT
		
		cudaDeviceSynchronize();
	}

	inline void TetMeshFEMShared::initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU)
	{
		TetMeshTopologyGPU* pTopology = getCPUTetMeshTopology(pTetMesh);
		pTetMeshGPU->pTopology = pTopology;

		pTetMeshGPU->tetRestVolume = pTetMesh->tetRestVolume.data();
		pTetMeshGPU->tetInvRestVolume = pTetMesh->tetInvRestVolume.data();
		pTetMeshGPU->vertexMass = pTetMesh->vertexMass.data();

		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetRestVolume);
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvRestVolume);
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertexInvMass);

		pTetMeshGPU->DmInvs = pTetMesh->DmInvs.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, DmInvs);

		pTetMeshGPU->vertPos = pTetMesh->mVertPos.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPos);
		pTetMeshGPU->vertPrevPos = pTetMesh->mVertPrevPos.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPrevPos);
		pTetMeshGPU->velocity = pTetMesh->mVelocity.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, velocity);
		//pTetMeshGPU->tetInvertedSign = pTetMesh->tetsInvertedSign.data();
		//pTetMeshGPU->verticesInvertedSign = pTetMesh->tetsInvertedSign.data();
		pTetMeshGPU->vertexFixedMask = pTetMesh->fixedMask.data();
	}
}
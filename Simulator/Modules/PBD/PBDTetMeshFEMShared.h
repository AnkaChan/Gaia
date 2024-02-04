#pragma once

#include "PBDTetMeshFEM.h"
#include "PBDTetMeshFEMGPU.h"
#include "CuMatrix/Buffers/ManagedBuffer.h"

#define SET_GPU_MESH_BUFER_PTR(gpuMesh, memberPtr)\
gpuMesh->memberPtr = memberPtr##Buffer->getGPUBuffer();

namespace GAIA {

	// a class manages and syncs data between CPU and GPU
	struct PBDTetMeshFEMShared {

		typedef std::shared_ptr<PBDTetMeshFEMShared> SharedPtr;
		typedef PBDTetMeshFEMShared* Ptr;

#ifdef USE_GPU_COLORING_DATA
		// those data can be passed by CPU, currently we don't need a copy of them on the GPU
		 ManagedBuffer<int32_t*>::SharedPtr tetsColoringCategoriesPointersBuffer;
		 ManagedBuffer<int32_t>::SharedPtr tetsColoringEachCategorySizeBuffer;
#endif // !USE_GPU_COLORING_DATA

		std::vector<ManagedBuffer<int32_t>::SharedPtr> tetsColoringEachCategoryBuffer;
		
		// managed buffers for non-variant data

		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetRestVolumeBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetInvRestVolumeBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertexInvMassBuffer;

		// flatten vertex data: 3 x nVerts
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPrevPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr velocityBuffer;

		// tet's vertex buffer: 4 x nTets
		ManagedBuffer<int32_t>::SharedPtr tetVIdsBuffer;
		// tet's inversion of restpose transformation: 9 x nTets
		ManagedBuffer<FloatingTypeGPU>::SharedPtr DmInvsBuffer;
		// tet's inversion sign buffer:  nTets
		ManagedBuffer<int8_t>::SharedPtr tetInvertedSignBuffer;
		// tet's inversion sign buffer:  nVerts
		ManagedBuffer<int8_t>::SharedPtr verticesInvertedSignBuffer;

#ifdef EVAL_TET_ACCESS_COUNT
		ManagedBuffer<int32_t>::SharedPtr tetAccessCountBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertAccessCountBuffer;
#endif //EVAL_TET_ACCESS_COUNT

		template <typename GPUTetMeshType>
		void initializeGPUTetMesh(PBDTetMeshFEM* pTetMesh, GPUTetMeshType* pTetMeshGPU);
		// Unlike the above one, this generates the GPU TetMesh with CPU buffers
		template <typename GPUTetMeshType>
		void initializeGPUTetMeshForDebug(PBDTetMeshFEM* pTetMesh, GPUTetMeshType* pTetMeshGPU);

		virtual void copyAllToGPU(bool sync = true, cudaStream_t stream = 0);

		virtual void syncToGPU(bool sync=true, cudaStream_t stream = 0);

		virtual void syncToCPU(bool sync = true, cudaStream_t stream = 0);

		void syncToCPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);
		void syncToGPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);

		virtual void syncToCPUInversionSignOnly(bool sync = true, cudaStream_t stream = 0);

		size_t numTetParallelizationGroups() {
			return tetsColoringEachCategoryBuffer.size();
		};

		int32_t* getTetParallelizationGroupBufferGPU(int32_t iGroup) {
			return tetsColoringEachCategoryBuffer[iGroup]->getGPUBuffer();
		};
		size_t getTetParallelizationGroupSize(int32_t iGroup) {
			return tetsColoringEachCategoryBuffer[iGroup]->getSize();
		};
	};

	template<typename GPUTetMeshType>
	inline void PBDTetMeshFEMShared::initializeGPUTetMesh(PBDTetMeshFEM* pTetMesh, GPUTetMeshType* pTetMeshGPU)
	{
		std::cout << "Initializing GPU tet mesh for base material.\n";
		tetRestVolumeBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->tetRestVolume.size(), true, pTetMesh->tetRestVolume.data());
		tetRestVolumeBuffer->toGPU(false);
		tetInvRestVolumeBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->tetInvRestVolume.size(), true, pTetMesh->tetInvRestVolume.data());
		tetInvRestVolumeBuffer->toGPU(false);
		vertexInvMassBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->vertexInvMass.size(), true, pTetMesh->vertexInvMass.data());
		vertexInvMassBuffer->toGPU(false);

		tetVIdsBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->tetVIds().size(), true, pTetMesh->tetVIds().data());
		tetVIdsBuffer->toGPU(false);

		DmInvsBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->DmInvs.size(), true, pTetMesh->DmInvs.data());
		DmInvsBuffer->toGPU(false);

		vertPosBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVertPos.size(), true, pTetMesh->mVertPos.data());
		vertPosBuffer->toGPU(false);
		vertPrevPosBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVertPrevPos.size(), true, pTetMesh->mVertPrevPos.data());
		vertPrevPosBuffer->toGPU(false);
		velocityBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->mVelocity.size(), true, pTetMesh->mVelocity.data());
		velocityBuffer->toGPU(false);

		tetInvertedSignBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->tetsInvertedSign.size(),
			true, pTetMesh->tetsInvertedSign.data());
		tetInvertedSignBuffer->toGPU(false);

		verticesInvertedSignBuffer = std::make_shared<ManagedBuffer<int8_t>>(pTetMesh->verticesInvertedSign.size(),
			true, pTetMesh->verticesInvertedSign.data());
		verticesInvertedSignBuffer->toGPU(false);

#ifdef USE_GPU_COLORING_DATA
		// intialize pointer buffer and each category buffer
		tetsColoringCategoriesPointersBuffer = std::make_shared<ManagedBuffer<int32_t*>>(pTetMesh->tetsColoringCategories().size(), true);
		tetsColoringEachCategorySizeBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->tetsColoringCategories().size(), true);
#endif // !USE_GPU_COLORING_DATA

		for (size_t iColorCategory = 0; iColorCategory < pTetMesh->tetsColoringCategories().size(); iColorCategory++)
		{
			tetsColoringEachCategoryBuffer.push_back(
				std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->tetsColoringCategories()[iColorCategory].size(), true, pTetMesh->tetsColoringCategories()[iColorCategory].data())
			);
			// copy CPU coloring buffer to GPU
			tetsColoringEachCategoryBuffer.back()->toGPU(false);

#ifdef USE_GPU_COLORING_DATA
			// set the size of this coloring category
			tetsColoringEachCategorySizeBuffer->getCPUBuffer()[iColorCategory] = tetsColoringEachCategoryBuffer.back()->getSize();
			// set the pointer to this coloring category's buffer
			tetsColoringCategoriesPointersBuffer->getCPUBuffer()[iColorCategory] = tetsColoringEachCategoryBuffer.back()->getGPUBuffer();
#endif // !USE_GPU_COLORING_DATA
		}

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

		pTetMeshGPU->nVerts = pTetMesh->numVertices();
		pTetMeshGPU->nTets = pTetMesh->numTets();

		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetRestVolume);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvRestVolume);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertexInvMass);

		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetVIds);

		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, DmInvs);

		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPos);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPrevPos);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, velocity);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvertedSign);
		SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, verticesInvertedSign);
			
		tetsColoringEachCategorySizeBuffer->toGPU();
		tetsColoringCategoriesPointersBuffer->toGPU();
		pTetMeshGPU->numTetsColoredCatergories = pTetMesh->tetsColoringCategories().size();
		pTetMeshGPU->tetsColoringEachCategorySize = tetsColoringEachCategorySizeBuffer->getGPUBuffer();
		pTetMeshGPU->tetsColoringCategoriesPointers = tetsColoringCategoriesPointersBuffer->getGPUBuffer();

#ifdef USE_GPU_COLORING_DATA
#endif // !USE_GPU_COLORING_DATA	

		pTetMeshGPU->dt = pTetMesh->dt;

		cudaDeviceSynchronize();
	}

	template<typename GPUTetMeshType>
	inline void PBDTetMeshFEMShared::initializeGPUTetMeshForDebug(PBDTetMeshFEM* pTetMesh, GPUTetMeshType* pTetMeshGPU)
	{
		pTetMeshGPU->nVerts = pTetMesh->numVertices();
		pTetMeshGPU->nTets = pTetMesh->numTets();

		pTetMeshGPU->tetRestVolume = pTetMesh->tetRestVolume.data();
		pTetMeshGPU->tetInvRestVolume = pTetMesh->tetInvRestVolume.data();
		pTetMeshGPU->vertexInvMass = pTetMesh->vertexInvMass.data();

		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetRestVolume);
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetInvRestVolume);
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertexInvMass);

		pTetMeshGPU->tetVIds = pTetMesh->tetVIds().data();
		// SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, tetVIds);

		pTetMeshGPU->DmInvs = pTetMesh->DmInvs.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, DmInvs);

		pTetMeshGPU->vertPos = pTetMesh->mVertPos.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPos);
		pTetMeshGPU->vertPrevPos = pTetMesh->mVertPrevPos.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, vertPrevPos);
		pTetMeshGPU->velocity = pTetMesh->mVelocity.data();
		//SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, velocity);
		pTetMeshGPU->tetInvertedSign = pTetMesh->tetsInvertedSign.data();
		pTetMeshGPU->verticesInvertedSign = pTetMesh->tetsInvertedSign.data();

		pTetMeshGPU->dt = pTetMesh->dt;
	}
}
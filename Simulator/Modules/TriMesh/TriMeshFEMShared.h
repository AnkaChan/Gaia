#pragma once

#include "TriMesh.h"

namespace GAIA {
	struct TriMeshTopologyShared {
		ManagedBuffer<int32_t>::SharedPtr faceVIdsBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborFacesBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborFaces_vertexOrderBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborFaces_infosBuffer;
		
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborBendingsBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborBendings_vertexOrderBuffer;
		ManagedBuffer<int32_t>::SharedPtr vertexNeighborBendings_infosBuffer;

		std::shared_ptr<DeviceClassBuffer<TriMeshTopologyGPU>> topologyGPUBuffer;
		std::shared_ptr<TriMeshTopologyGPU> topologyGPUForCPUDebug;

		void initialize(TriMeshTopology* pTopology);
		void initializeForCPUDebug(TriMeshTopology* pTopology);

	};


	struct TriMeshFEMShared
	{
		typedef std::shared_ptr<TriMeshFEMShared> SharedPtr;
		typedef TriMeshFEMShared* Ptr;

		// managed buffers for non-variant data
		ManagedBuffer<FloatingTypeGPU>::SharedPtr faceRestAreaBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr faceInvRestAreaBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertexMassBuffer;

		// data arranged per face
		ManagedBuffer<FloatingTypeGPU>::SharedPtr DmInvs;            // tet's mats, 9 x nTets
		ManagedBuffer<FloatingTypeGPU>::SharedPtr faceRestArea;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr faceInvRestArea;

		// flatten vertex data: 3 x nVerts
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr vertPrevPosBuffer;
		ManagedBuffer<FloatingTypeGPU>::SharedPtr velocityBuffer;

		// from the mesh path to trimesh GPU topologies on GPU memory
		static std::map<std::string, std::shared_ptr<TriMeshTopologyShared>> topologies;
		static std::mutex topologies_lock;
		TriMeshTopologyGPU* getGPUTetMeshTopology(TriMeshFEM* pTetMesh);
		TriMeshTopologyGPU* getCPUTetMeshTopology(TriMeshFEM* pTetMesh);

		// Generate GPU buffers, then copy them to a GPU mesh on CPU, which is ready to be copied to GPU
		void initializeGPUTetMesh(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);
		// Unlike the above one, this generates the GPU TetMesh with CPU buffers
		void initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU);


		virtual void copyAllToGPU(bool sync = true, cudaStream_t stream = 0);
		// sync all the mutable arrays, namely vertPos, vertPrevPos, and velocity to GPU
		virtual void syncToGPU(bool sync = true, cudaStream_t stream = 0);
		// sync all the mutable arrays, namely vertPos, vertPrevPos, and velocity to CPU
		virtual void syncToCPU(bool sync = true, cudaStream_t stream = 0);

		void syncToCPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);
		void syncToGPUVertPosOnly(bool sync = true, cudaStream_t stream = 0);

		virtual void syncToCPUInversionSignOnly(bool sync = true, cudaStream_t stream = 0);
	};
}
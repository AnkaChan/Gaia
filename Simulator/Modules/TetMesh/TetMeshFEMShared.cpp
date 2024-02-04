#include "TetMeshFEMShared.h"

using namespace GAIA;

//void GAIA::TetMeshFEMShared::createGPUTetMesh()
//{
//	std::cout << "Creating GPU tet mesh for the base material.\n";
//	pTetMeshGPU = std::make_shared<ManagedClassBuffer<TetMeshFEMGPU>>();
//}
std::map<std::string, std::shared_ptr<TetMeshTopologyShared>> TetMeshFEMShared::topologies;
std::mutex TetMeshFEMShared::topologies_lock;


TetMeshTopologyGPU* GAIA::TetMeshFEMShared::getGPUTetMeshTopology(TetMeshFEM* pTetMesh)
{
	bool alreadyComputed = false;
	{
		std::lock_guard<std::mutex> topologyLockGuard(topologies_lock);
		std::string modelPath = pTetMesh->pObjectParams->path;


		auto pTopoItem = topologies.find(modelPath);

		if (pTopoItem == topologies.end())
		{
			// make a GPU topology and copy to the GPU memory
			auto pTopologyShared = std::make_shared<TetMeshTopologyShared>();
			pTopologyShared->initialize(pTetMesh->pTopology.get());
			topologies.insert({ modelPath, pTopologyShared });
			return pTopologyShared->topologyGPUBuffer->getData();
		}
		else
		{
			TetMeshTopologyGPU* pTopologyGPU = pTopoItem->second->topologyGPUBuffer->getData();
			return pTopologyGPU;
		}
	}
}

TetMeshTopologyGPU* GAIA::TetMeshFEMShared::getCPUTetMeshTopology(TetMeshFEM* pTetMesh)
{
	bool alreadyComputed = false;
	{
		std::lock_guard<std::mutex> topologyLockGuard(topologies_lock);
		std::string modelPath = pTetMesh->pObjectParams->path;


		auto pTopoItem = topologies.find(modelPath);

		if (pTopoItem == topologies.end())
		{
			// make a GPU topology and copy to the GPU memory
			auto pTopologyShared = std::make_shared<TetMeshTopologyShared>();
			pTopologyShared->initializeForCPUDebug(pTetMesh->pTopology.get());
			topologies.insert({ modelPath, pTopologyShared });
			return pTopologyShared->topologyGPUForCPUDebug.get();
		}
		else
		{
			if (pTopoItem->second->topologyGPUForCPUDebug.get() == nullptr)
			{
				pTopoItem->second->initializeForCPUDebug(pTetMesh->pTopology.get());
			}
			TetMeshTopologyGPU* pTopologyGPU = pTopoItem->second->topologyGPUForCPUDebug.get();
			return pTopologyGPU;
		}
	}
}


void GAIA::TetMeshFEMShared::copyAllToGPU(bool sync, cudaStream_t stream)
{
	tetRestVolumeBuffer->toGPU(false);
	tetInvRestVolumeBuffer->toGPU(false);
	vertexMassBuffer->toGPU(false);

	DmInvsBuffer->toGPU(false);

	vertPosBuffer->toGPU(false);
	vertPrevPosBuffer->toGPU(false);
	velocityBuffer->toGPU(false);

	//tetInvertedSignBuffer->toGPU(false);
	//verticesInvertedSignBuffer->toGPU(false);

	cudaDeviceSynchronize();

}

void GAIA::TetMeshFEMShared::syncToGPU(bool sync , cudaStream_t stream)
{
	vertPosBuffer->toGPU(false, stream);
	vertPrevPosBuffer->toGPU(false, stream);
	velocityBuffer->toGPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::TetMeshFEMShared::syncToCPU(bool sync, cudaStream_t stream)
{
	vertPosBuffer->toCPU(false, stream);
	vertPrevPosBuffer->toCPU(false, stream);
	velocityBuffer->toCPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::TetMeshFEMShared::syncToCPUVertPosOnly(bool sync, cudaStream_t stream)
{
	vertPosBuffer->toCPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::TetMeshFEMShared::syncToGPUVertPosOnly(bool sync, cudaStream_t stream)
{
	vertPosBuffer->toGPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::TetMeshFEMShared::syncToCPUInversionSignOnly(bool sync, cudaStream_t stream)
{
	//tetInvertedSignBuffer->toCPU(false, stream);
	//verticesInvertedSignBuffer->toCPU(false, stream);

	if (sync)
	{
		cudaStreamSynchronize(stream);
	}
}

void GAIA::TetMeshTopologyShared::initialize(TetMeshTopology* pTopology)
{
	TetMeshTopologyGPU topologyOnCPU;
	TetMeshTopologyGPU* pTopologyOnCPU = &topologyOnCPU;
	
	topologyOnCPU.nTets = pTopology->numTets();
	topologyOnCPU.nVerts = pTopology->numVertices();

	tetVIdsBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTopology->tetVIds.size(), true, pTopology->tetVIds.data());
	tetVIdsBuffer->toGPU(false);
	SET_GPU_MESH_BUFER_PTR(pTopologyOnCPU, tetVIds);

	vertexNeighborTetsBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTopology->vertexNeighborTets.size(), true, pTopology->vertexNeighborTets.data());
	vertexNeighborTetsBuffer->toGPU(false);
	SET_GPU_MESH_BUFER_PTR(pTopologyOnCPU, vertexNeighborTets);

	vertexNeighborTets_vertexOrderBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTopology->vertexNeighborTets_vertexOrder.size(), true, pTopology->vertexNeighborTets_vertexOrder.data());
	vertexNeighborTets_vertexOrderBuffer->toGPU(false);
	SET_GPU_MESH_BUFER_PTR(pTopologyOnCPU, vertexNeighborTets_vertexOrder);

	vertexNeighborTets_infosBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTopology->vertexNeighborTets_infos.size(), true, pTopology->vertexNeighborTets_infos.data());
	vertexNeighborTets_infosBuffer->toGPU(false);
	SET_GPU_MESH_BUFER_PTR(pTopologyOnCPU, vertexNeighborTets_infos);

	topologyGPUBuffer = std::make_shared<DeviceClassBuffer<TetMeshTopologyGPU>>();
	topologyGPUBuffer->fromCPU(pTopologyOnCPU);
}

void GAIA::TetMeshTopologyShared::initializeForCPUDebug(TetMeshTopology* pTopology)
{
	topologyGPUForCPUDebug = std::make_shared<TetMeshTopologyGPU>();

	topologyGPUForCPUDebug->nTets = pTopology->numTets();
	topologyGPUForCPUDebug->nVerts = pTopology->numVertices();

	topologyGPUForCPUDebug->tetVIds = pTopology->tetVIds.data();
	topologyGPUForCPUDebug->vertexNeighborTets = pTopology->vertexNeighborTets.data();
	topologyGPUForCPUDebug->vertexNeighborTets_vertexOrder = pTopology->vertexNeighborTets_vertexOrder.data();
	topologyGPUForCPUDebug->vertexNeighborTets_infos = pTopology->vertexNeighborTets_infos.data();




}

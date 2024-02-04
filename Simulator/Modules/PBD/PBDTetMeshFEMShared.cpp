#include "PBDTetMeshFEMShared.h"



//void GAIA::TetMeshFEMShared::createGPUTetMesh()
//{
//	std::cout << "Creating GPU tet mesh for the base material.\n";
//	pTetMeshGPU = std::make_shared<ManagedClassBuffer<TetMeshFEMGPU>>();
//}


void GAIA::PBDTetMeshFEMShared::copyAllToGPU(bool sync, cudaStream_t stream)
{
	tetRestVolumeBuffer->toGPU(false);
	tetInvRestVolumeBuffer->toGPU(false);
	vertexInvMassBuffer->toGPU(false);

	tetVIdsBuffer->toGPU(false);
	DmInvsBuffer->toGPU(false);

	vertPosBuffer->toGPU(false);
	vertPrevPosBuffer->toGPU(false);
	velocityBuffer->toGPU(false);

	tetInvertedSignBuffer->toGPU(false);
	verticesInvertedSignBuffer->toGPU(false);

	for (size_t iColorCategory = 0; iColorCategory < tetsColoringEachCategoryBuffer.size(); iColorCategory++)
	{
		tetsColoringEachCategoryBuffer[iColorCategory]->toGPU(false);
	}

	tetsColoringEachCategorySizeBuffer->toGPU(false);
	tetsColoringCategoriesPointersBuffer->toGPU(false);

	cudaDeviceSynchronize();

}

void GAIA::PBDTetMeshFEMShared::syncToGPU(bool sync , cudaStream_t stream)
{
	vertPosBuffer->toGPU(false, stream);
	vertPrevPosBuffer->toGPU(false, stream);
	velocityBuffer->toGPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::PBDTetMeshFEMShared::syncToCPU(bool sync, cudaStream_t stream)
{


	vertPosBuffer->toCPU(false, stream);
	vertPrevPosBuffer->toCPU(false, stream);
	velocityBuffer->toCPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::PBDTetMeshFEMShared::syncToCPUVertPosOnly(bool sync, cudaStream_t stream)
{
	vertPosBuffer->toCPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::PBDTetMeshFEMShared::syncToGPUVertPosOnly(bool sync, cudaStream_t stream)
{
	vertPosBuffer->toGPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}
}

void GAIA::PBDTetMeshFEMShared::syncToCPUInversionSignOnly(bool sync, cudaStream_t stream)
{
	tetInvertedSignBuffer->toCPU(false, stream);
	verticesInvertedSignBuffer->toCPU(false, stream);

	if (sync)
	{
		cudaStreamSynchronize(stream);
	}
}


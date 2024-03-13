#include "VBDPhysics.h"
#include "LineSearchUtilities.h"

void GAIA::LineSearchUtilities::initialize(const VBDPhysics& physics)
{
	const std::vector<VBDBaseTetMesh::SharedPtr>& tMeshes = physics.tMeshes;
	size_t numAllTets = physics.numAllTets;

	tetElasticEnergyBuffer = std::make_shared<ManagedBuffer<float>>(numAllTets, true);
	tetAllParallelGroupsBuffer = std::make_shared<ManagedBuffer<int32_t>>(2 * numAllTets, true);
	size_t tetAccumulate = 0;
	for (int32_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		const size_t numTets = pTetMesh->numTets();
		for (int32_t iTet = 0; iTet < numTets; iTet++)
		{
			tetAllParallelGroupsBuffer->getCPUBuffer()[2 * (iTet + tetAccumulate)] = iMesh;
			tetAllParallelGroupsBuffer->getCPUBuffer()[2 * (iTet + tetAccumulate) + 1] = iTet;
		}
		tetAccumulate += numTets;
	}
	assert(tetAccumulate == numAllTets);
	tetAllParallelGroupsBuffer->toGPU();
}

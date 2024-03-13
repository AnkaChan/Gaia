#include "VBDPhysics.h"
#include "GD_SolverUtilities.h"

void GAIA::GD_SolverUtilities::initialize(VBDPhysics& physics)
{
	auto& tMeshes = physics.tMeshes;
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		previousPositionsForLineSearchBuffer.push_back(
			std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->vertices().size(), false)
		);
		prevPrevPositionsForAcceleratorBuffer.push_back(
			std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->vertices().size(), true)
		);
	}

	stepSizePrevStep = physics.physicsParams().stepSizeGD;
}

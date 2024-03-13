#pragma once
#include "../Types/Types.h"
#include "VBD_BaseMaterial.h"
namespace GAIA {
	struct VBDPhysics;

	struct GD_SolverUtilities {
		typedef std::shared_ptr<GD_SolverUtilities> SharedPtr;
		typedef GD_SolverUtilities* Ptr;

		void initialize(VBDPhysics& physics);
		std::vector<ManagedBuffer<FloatingTypeGPU>::SharedPtr> previousPositionsForLineSearchBuffer;
		std::vector<ManagedBuffer<FloatingTypeGPU>::SharedPtr> prevPrevPositionsForAcceleratorBuffer;

		FloatingType stepSizePrevStep = -1.f;
		FloatingType omegaPrev = 1.0f;
	};
}
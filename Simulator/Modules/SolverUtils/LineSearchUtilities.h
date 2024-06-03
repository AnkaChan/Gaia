#pragma once
#include "../Types/Types.h"
#include "../VBD/VBD_BaseMaterial.h"

namespace GAIA {
	struct VBDPhysics;
	struct LineSearchUtilities {
		typedef std::shared_ptr<LineSearchUtilities> SharedPtr;
		typedef LineSearchUtilities* Ptr;

		void initialize(const VBDPhysics& physics);

		void recordPrevEnergy() {
			ePrev = e;
			eInertiaPrev = eInertia;
			eElasticPrev = eElastic;
		}

		ManagedBuffer<FloatingTypeGPU>::SharedPtr tetElasticEnergyBuffer;
		ManagedBuffer<int32_t>::SharedPtr tetAllParallelGroupsBuffer;


		FloatingType ePrev = 0.f;
		FloatingType eInertiaPrev = 0.f;
		FloatingType eElasticPrev = 0.f;

		FloatingType e = 0.f;
		FloatingType eInertia = 0.f;
		FloatingType eElastic = 0.f;
	};

}
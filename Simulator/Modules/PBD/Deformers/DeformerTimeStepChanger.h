#pragma once
#include <functional>
#include <memory>

#include "../Deformer.h"
#include "../PBDPhysics.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

namespace EBD {
	struct DeformerTimeStepChanger : public BaseDeformer {
		DeformerTimeStepChanger(PBDPhysics& physics, int targetNumSubSteps, int targetNumIterations,
			FloatingType in_deformationEndTime, FloatingType in_deformationStartTime = 0)
			: numSubSteps1(targetNumSubSteps)
			, numIteration1(targetNumIterations)
			, deformationEndTime(in_deformationEndTime)
			, deformationStartTime(in_deformationStartTime)
		{
			numSubSteps0 = physics.physicsParams().numSubsteps;
			numIteration0 = physics.physicsParams().iterations;
		}

		virtual void operator()(PBDPhysics* pPBDPhysics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
			if (curTime < deformationStartTime)
				return;
			FloatingType w = (curTime - deformationStartTime) / (deformationEndTime - deformationStartTime);

			if (w > 1.0) w = 1.0;

			int newNumSubSteps = int(std::round(numSubSteps0 * (1 - w) + numSubSteps1 * w));
			int newIterations = int(std::round(numIteration0 * (1 - w) + numIteration1 * w));
			if (pPBDPhysics->physicsParams().numSubsteps != newNumSubSteps) {
				pPBDPhysics->physicsParams().numSubsteps = newNumSubSteps;
				std::cout << "Num subteps changed to: " << pPBDPhysics->physicsParams().numSubsteps << "\n";

				pPBDPhysics->changeNumSubSteps(newNumSubSteps);
			}

			if (pPBDPhysics->physicsParams().iterations != newIterations) {
				pPBDPhysics->physicsParams().iterations = newIterations;
				std::cout << "Num iterations changed to: " << pPBDPhysics->physicsParams().iterations << "\n";

			}
		};

		
		int numSubSteps0;
		int numSubSteps1;

		int numIteration0;
		int numIteration1;
		FloatingType deformationEndTime;
		FloatingType deformationStartTime;
	};


}
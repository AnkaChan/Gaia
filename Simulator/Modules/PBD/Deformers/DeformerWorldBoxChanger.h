#pragma once
#include <functional>
#include <memory>

#include "../Deformer.h"
#include "../PBDPhysics.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

namespace EBD {
	struct DeformerWorldBoxChanger : public BaseDeformer {
		DeformerWorldBoxChanger(const Eigen::Matrix<FloatingType, 3, 2>& in_targetWorldBox, const Eigen::Matrix<FloatingType, 3, 2>& in_originalWorldBox,
			FloatingType in_deformationEndTime, FloatingType in_deformationStartTime = 0)
			: targetWorldBox(in_targetWorldBox)
			, originalWorldBox(in_originalWorldBox)
			, deformationStartTime(in_deformationStartTime)
			, deformationEndTime(in_deformationEndTime)
		{


		}

		virtual void operator()(PBDPhysics* pPBDPhysics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {

			FloatingType w = (curTime - deformationStartTime) / (deformationEndTime - deformationStartTime);

			if (w > 1.0) w = 1.0;

			pPBDPhysics->physicsParams().worldBounds = originalWorldBox * (1 - w) + targetWorldBox * w;
			if (iSubstep == 0 && iIter == 0) {
				std::cout << "World box: \n" << pPBDPhysics->physicsParams().worldBounds.transpose() << "\n";
			}
		};

		Eigen::Matrix<FloatingType, 3, 2> targetWorldBox;
		Eigen::Matrix<FloatingType, 3, 2> originalWorldBox;
		FloatingType deformationEndTime;
		FloatingType deformationStartTime;

	};

	struct DeformerWorldBoxChanger2Stages : public BaseDeformer {
		DeformerWorldBoxChanger2Stages(const Eigen::Matrix<FloatingType, 3, 2>& in_targetWorldBox, const Eigen::Matrix<FloatingType, 3, 2>& in_originalWorldBox,
			FloatingType in_linearDeformationStartTime, // FloatingType in_linearDeformationEndTime, // linear deformation ends when exponential deformation ends
			FloatingType in_linearDeformationEndRatio,  // blending ratio when the linear deformation stops
			FloatingType in_exponentialDeformationStartTime, FloatingType in_exponentialDeformationEndTime, 
			FloatingType in_decreaseRationPerSec)
			: targetWorldBox(in_targetWorldBox)
			, originalWorldBox(in_originalWorldBox)
			, linearDeformationStartTime(in_linearDeformationStartTime)
			, linearDeformationEndRatio(in_linearDeformationEndRatio)
			, exponentialDeformationStartTime(in_exponentialDeformationStartTime)
			, exponentialDeformationEndTime(in_exponentialDeformationEndTime)
			, decreaseRatioPerSec(in_decreaseRationPerSec)
		{

		}

		virtual void operator()(PBDPhysics* pPBDPhysics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
			if (iIter != 0)
			{
				return;
			}
			if (curTime < exponentialDeformationStartTime)
			{
				FloatingType w = (curTime - linearDeformationStartTime) / (exponentialDeformationStartTime - linearDeformationStartTime);

				if (w > 1.0) w = 1.0;

				w = linearDeformationEndRatio * w;
				pPBDPhysics->physicsParams().worldBounds = originalWorldBox * (1 - w) + targetWorldBox * w;
				if (iSubstep == 0 && iIter == 0) {
					std::cout << "Linear deformation stage. ";
				}
			}
			else {

				FloatingType w = (1.f - linearDeformationEndRatio)  * pow(decreaseRatioPerSec, curTime - exponentialDeformationStartTime);
				pPBDPhysics->physicsParams().worldBounds = originalWorldBox * w + targetWorldBox * (1 - w);
				if (iSubstep == 0 && iIter == 0) {
					std::cout << "Exponential deformation stage. ";
				}
			}

			if (iSubstep == 0 && iIter == 0) {
				std::cout << "World box: \n" << pPBDPhysics->physicsParams().worldBounds.transpose() << "\n";
			}
		};

		Eigen::Matrix<FloatingType, 3, 2> targetWorldBox;
		Eigen::Matrix<FloatingType, 3, 2> originalWorldBox;
		FloatingType linearDeformationStartTime;
		FloatingType linearDeformationEndRatio;

		FloatingType exponentialDeformationEndTime;
		FloatingType exponentialDeformationStartTime;
		FloatingType decreaseRatioPerSec;

	};
}
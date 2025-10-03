#pragma once

#include <Framework/BaseClothSimPhsicsFramework.h>

namespace GAIA {

	struct VBDClothPhysicsParameters : public BasePhysicsParams {
		typedef std::shared_ptr<VBDClothPhysicsParameters> SharedPtr;
		typedef VBDClothPhysicsParameters* Ptr;

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);

		bool fit = true;
		int maxNumIteration = 10000;
		int evaluationSteps = 1;
		FloatingType stepSize = 1.0f;
		FloatingType stepSizeGD = 1e-3f;
		FloatingType minStepSizeGD = 1e-3f;
		bool applyAcceleration = false;
		FloatingType acceleratorRho = 0.95f;

		// for convergence evaluation
		bool evaluateConvergence = true;
		FloatingType convergenceEnergyChangeThres = 1e-6f;
		FloatingType convergenceAvgNormThres = 20.f;
		FloatingType convergenceAvgNormChangeThres = 1e-4f;

		bool useParallel = true;

		// 1: avgGradNorm + avgGradNormChange, 1: energy
		int convergenceType = 0;

		// initialization
		bool associateGravityWithInertia = true;
		int initializationType = 0; // 0: with inertia, 1: previous location (to ensure the initial state is collision free)
		// if false, will set the gradient wrt. intertia to zero, in this case it will become a quasi-static fit
		bool withInertia = true;

		// optimization method type:
		bool useNewton = false;
		int newtonSolverType = 0; // 0: direct, 1: CG

		// for back tracing line search
		bool useLineSearch = false;
		FloatingType c = 0.01f;
		FloatingType tau = 0.5f;
		int backtracingLineSearchMaxIters = 8;
		int numLineSearchSteps = 2;
		int lineSearchEveryNumSteps = -1;

		std::vector<std::string> targetFiles;


		// closest point query
		int closestPointQuerySteps = 1;
		FloatingType registrationEnergyWeight = 1e8f;

		// collision detecion
		bool handleCollision = true;
		FloatingType contactStiffness = 1e4f;
		FloatingType contactRadius = 1.f;
		FloatingType thickness = 0.0f;
		FloatingType conservativeStepRelaxation = 0.45f;

		FloatingType boundaryContactStiffness = 1e4f;
		FloatingType boundaryFrictionEpsValue = 1e-2f;

		int contactDetectionIters = 1;
		int contactEnergyTypes = 1; // 0: quadratic, 1: quadratic-logrithmic-2-stages
		int contactBVHReconstructionSteps = 32;

		// material solve
		FloatingType degenerateTriangleThres = 1e-6f;

		// for debug
		bool saveIntermediateResults = false;
		int debugOutputInterval = 1;


	};

	inline bool VBDClothPhysicsParameters::fromJson(nlohmann::json& objectParam)
	{
		BasePhysicsParams::fromJson(objectParam);
		EXTRACT_FROM_JSON(objectParam, fit);
		EXTRACT_FROM_JSON(objectParam, useParallel);
		EXTRACT_FROM_JSON(objectParam, maxNumIteration);
		EXTRACT_FROM_JSON(objectParam, evaluationSteps);
		EXTRACT_FROM_JSON(objectParam, stepSize);
		EXTRACT_FROM_JSON(objectParam, minStepSizeGD);
		EXTRACT_FROM_JSON(objectParam, associateGravityWithInertia);
		EXTRACT_FROM_JSON(objectParam, useNewton);
		EXTRACT_FROM_JSON(objectParam, newtonSolverType);
		EXTRACT_FROM_JSON(objectParam, acceleratorRho);
		EXTRACT_FROM_JSON(objectParam, applyAcceleration);

		EXTRACT_FROM_JSON(objectParam, useLineSearch);
		EXTRACT_FROM_JSON(objectParam, c);
		EXTRACT_FROM_JSON(objectParam, tau);
		EXTRACT_FROM_JSON(objectParam, backtracingLineSearchMaxIters);
		EXTRACT_FROM_JSON(objectParam, numLineSearchSteps);
		EXTRACT_FROM_JSON(objectParam, lineSearchEveryNumSteps);


		EXTRACT_FROM_JSON(objectParam, evaluateConvergence);
		EXTRACT_FROM_JSON(objectParam, convergenceType);
		EXTRACT_FROM_JSON(objectParam, convergenceEnergyChangeThres);
		EXTRACT_FROM_JSON(objectParam, convergenceAvgNormThres);
		EXTRACT_FROM_JSON(objectParam, convergenceAvgNormChangeThres);

		EXTRACT_FROM_JSON(objectParam, contactStiffness);
		EXTRACT_FROM_JSON(objectParam, contactRadius);
		EXTRACT_FROM_JSON(objectParam, contactBVHReconstructionSteps);
		EXTRACT_FROM_JSON(objectParam, conservativeStepRelaxation);
		EXTRACT_FROM_JSON(objectParam, thickness);
		EXTRACT_FROM_JSON(objectParam, handleCollision);
		EXTRACT_FROM_JSON(objectParam, boundaryFrictionEpsValue);

		EXTRACT_FROM_JSON(objectParam, targetFiles);
		EXTRACT_FROM_JSON(objectParam, withInertia);
		EXTRACT_FROM_JSON(objectParam, registrationEnergyWeight);

		EXTRACT_FROM_JSON(objectParam, handleCollision);

		EXTRACT_FROM_JSON(objectParam, saveIntermediateResults);
		EXTRACT_FROM_JSON(objectParam, debugOutputInterval);

		return true;
	}

	inline bool VBDClothPhysicsParameters::toJson(nlohmann::json& objectParam)
	{
		BasePhysicsParams::toJson(objectParam);
		PUT_TO_JSON(objectParam, fit);
		PUT_TO_JSON(objectParam, useParallel);
		PUT_TO_JSON(objectParam, maxNumIteration);
		PUT_TO_JSON(objectParam, evaluationSteps);
		PUT_TO_JSON(objectParam, stepSize);
		PUT_TO_JSON(objectParam, minStepSizeGD);
		PUT_TO_JSON(objectParam, associateGravityWithInertia);
		PUT_TO_JSON(objectParam, useNewton);
		PUT_TO_JSON(objectParam, newtonSolverType);
		PUT_TO_JSON(objectParam, acceleratorRho);
		PUT_TO_JSON(objectParam, applyAcceleration);

		PUT_TO_JSON(objectParam, useLineSearch);
		PUT_TO_JSON(objectParam, c);
		PUT_TO_JSON(objectParam, tau);
		PUT_TO_JSON(objectParam, backtracingLineSearchMaxIters);
		PUT_TO_JSON(objectParam, numLineSearchSteps);
		PUT_TO_JSON(objectParam, lineSearchEveryNumSteps);

		PUT_TO_JSON(objectParam, evaluateConvergence);
		PUT_TO_JSON(objectParam, convergenceType);
		PUT_TO_JSON(objectParam, convergenceEnergyChangeThres);
		PUT_TO_JSON(objectParam, convergenceAvgNormThres);
		PUT_TO_JSON(objectParam, convergenceAvgNormChangeThres);

		PUT_TO_JSON(objectParam, contactStiffness);
		PUT_TO_JSON(objectParam, contactRadius);
		PUT_TO_JSON(objectParam, contactBVHReconstructionSteps);
		PUT_TO_JSON(objectParam, conservativeStepRelaxation);
		PUT_TO_JSON(objectParam, thickness);
		PUT_TO_JSON(objectParam, handleCollision);
		PUT_TO_JSON(objectParam, boundaryFrictionEpsValue);

		PUT_TO_JSON(objectParam, targetFiles);
		PUT_TO_JSON(objectParam, withInertia);
		PUT_TO_JSON(objectParam, registrationEnergyWeight);

		PUT_TO_JSON(objectParam, handleCollision);

		PUT_TO_JSON(objectParam, saveIntermediateResults);
		PUT_TO_JSON(objectParam, debugOutputInterval);

		return true;
	}
}
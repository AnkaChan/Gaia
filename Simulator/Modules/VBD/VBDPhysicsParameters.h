#pragma once
#include "../Parameters/PhysicsParameters.h"

namespace GAIA {
	struct VBDPhysicsParameters : public BasePhysicsParams {
		typedef std::shared_ptr<VBDPhysicsParameters> SharedPtr;
		typedef VBDPhysicsParameters* Ptr;

		FloatingType stepSize = 1.0f;
		FloatingType stepSizeGD = 1e-3f;
		FloatingType convergenceEnergyChangeThres = 1e-6f;
		FloatingType convergenceAvgNormThres = 1e-2f;
		FloatingType convergenceAvgNormChangeThres = 1e-4f;

		// 1: avgGradNorm + avgGradNormChange, 1: energy
		int convergenceType = 0;

		// initialization
		bool associateGravityWithInertia = true;
		int initializationType = 0; // 0: with inertia, 1: previous location (to ensure the initial state is collision free)

		// for back tracing line search
		FloatingType c = 0.01f;
		FloatingType tau = 0.5f;

		FloatingType c_GD = 0.01f;
		FloatingType tau_GD = 0.1f;

		int backtracingLineSearchMaxIters = 16;
		FloatingType backtracingLineSearchAlpha = 1.0; // initial step size
		FloatingType backtracingLineSearchTau = 0.7; // step size multiplier
		FloatingType backtracingLineSearchC = 0.0; // first wolfe condition multiplier

		// collision detecion
		bool handleCollision = false;
		FloatingType collisionStiffness = 1e5f;
		FloatingType collisionAirDistance = 0.0f;
		int contactDetectionIters = 1;
		int collisionEnergyType = 1; // 0: point-point, 1: point-plane
		int collisionSolutionType = 0; // 0: serial, 1: hybrid
		int intermediateCollisionIterations = -1;
		FloatingType activeCollisionListPreAllocationRatio = 1.0f;

		FloatingType boundaryCollisionStiffness = 1e5f;
		FloatingType boundaryFrictionEpsV = 1.0f;
		int frictionStartIter = 0;

		// material solve
		FloatingType degenerateTriangleThres = 1e-6f;

		// solver
		bool useGPU = false;
		bool useDouble3x3 = false; // only for CPU
		bool useNewton = false;
		bool useGDSolver = false;  // only for GPU
		bool GDSolverUseBlockJacobi = false;  // only for GPU
		bool useLineSearch = false;

		bool NewtonUseCG = false;

		int lineSearchGapIter = 8;

		// accelerator
		bool useAccelerator = false;
		FloatingType acceleratorRho = 0.93f;

		// for debug
		bool saveIntermediateResults = false;
		int debugOutputInterval = 1;
		FloatingType collisionOffHeight = 1e7;
		bool evaluateConvergence = false;
		bool saveConvergenceEvaluationResults = false;
		int evaluationSteps = 10;

		// run GPU code on CPU
		bool debugGPU = false;

		// parallelism
		int numThreadsVBDSolve = 16;

		virtual bool fromJson(nlohmann::json& physicsParams);
		virtual bool toJson(nlohmann::json& physicsParams);
	};

	inline bool GAIA::VBDPhysicsParameters::fromJson(nlohmann::json& physicsParams)
	{
		BasePhysicsParams::fromJson(physicsParams);

		EXTRACT_FROM_JSON(physicsParams, stepSize);
		EXTRACT_FROM_JSON(physicsParams, stepSizeGD);

		EXTRACT_FROM_JSON(physicsParams, convergenceEnergyChangeThres);
		EXTRACT_FROM_JSON(physicsParams, convergenceAvgNormThres);
		EXTRACT_FROM_JSON(physicsParams, convergenceAvgNormChangeThres);

		EXTRACT_FROM_JSON(physicsParams, convergenceType);

		EXTRACT_FROM_JSON(physicsParams, associateGravityWithInertia);
		EXTRACT_FROM_JSON(physicsParams, initializationType);

		EXTRACT_FROM_JSON(physicsParams, c);
		EXTRACT_FROM_JSON(physicsParams, tau);

		EXTRACT_FROM_JSON(physicsParams, c_GD);
		EXTRACT_FROM_JSON(physicsParams, tau_GD);

		EXTRACT_FROM_JSON(physicsParams, backtracingLineSearchMaxIters);
		EXTRACT_FROM_JSON(physicsParams, backtracingLineSearchAlpha);
		EXTRACT_FROM_JSON(physicsParams, backtracingLineSearchTau);
		EXTRACT_FROM_JSON(physicsParams, backtracingLineSearchC);

		EXTRACT_FROM_JSON(physicsParams, useDouble3x3);
		EXTRACT_FROM_JSON(physicsParams, frictionStartIter);

		EXTRACT_FROM_JSON(physicsParams, useGPU);
		EXTRACT_FROM_JSON(physicsParams, useGDSolver);
		EXTRACT_FROM_JSON(physicsParams, useNewton);
		EXTRACT_FROM_JSON(physicsParams, useLineSearch);
		EXTRACT_FROM_JSON(physicsParams, lineSearchGapIter);
		EXTRACT_FROM_JSON(physicsParams, useAccelerator);
		EXTRACT_FROM_JSON(physicsParams, acceleratorRho);
		EXTRACT_FROM_JSON(physicsParams, GDSolverUseBlockJacobi);
		EXTRACT_FROM_JSON(physicsParams, NewtonUseCG);

		// collision detecion
		EXTRACT_FROM_JSON(physicsParams, handleCollision);
		EXTRACT_FROM_JSON(physicsParams, collisionStiffness);
		EXTRACT_FROM_JSON(physicsParams, contactDetectionIters);
		EXTRACT_FROM_JSON(physicsParams, collisionEnergyType);
		EXTRACT_FROM_JSON(physicsParams, collisionSolutionType);
		EXTRACT_FROM_JSON(physicsParams, collisionAirDistance);
		EXTRACT_FROM_JSON(physicsParams, boundaryCollisionStiffness);
		EXTRACT_FROM_JSON(physicsParams, intermediateCollisionIterations);
		EXTRACT_FROM_JSON(physicsParams, boundaryFrictionEpsV);
		EXTRACT_FROM_JSON(physicsParams, activeCollisionListPreAllocationRatio);
		// for debug
		EXTRACT_FROM_JSON(physicsParams, saveIntermediateResults);
		EXTRACT_FROM_JSON(physicsParams, debugOutputInterval);
		EXTRACT_FROM_JSON(physicsParams, collisionOffHeight);
		EXTRACT_FROM_JSON(physicsParams, evaluateConvergence);
		EXTRACT_FROM_JSON(physicsParams, evaluationSteps);
		EXTRACT_FROM_JSON(physicsParams, saveConvergenceEvaluationResults);

		EXTRACT_FROM_JSON(physicsParams, numThreadsVBDSolve);

		return true;
	}
	inline bool VBDPhysicsParameters::toJson(nlohmann::json& physicsParams)
	{
		BasePhysicsParams::toJson(physicsParams);

		PUT_TO_JSON(physicsParams, evaluationSteps);
		PUT_TO_JSON(physicsParams, stepSize);
		PUT_TO_JSON(physicsParams, stepSizeGD);

		PUT_TO_JSON(physicsParams, convergenceEnergyChangeThres);
		PUT_TO_JSON(physicsParams, convergenceAvgNormThres);
		PUT_TO_JSON(physicsParams, convergenceAvgNormChangeThres);

		PUT_TO_JSON(physicsParams, convergenceType);

		PUT_TO_JSON(physicsParams, associateGravityWithInertia);
		PUT_TO_JSON(physicsParams, initializationType);

		PUT_TO_JSON(physicsParams, c);
		PUT_TO_JSON(physicsParams, tau);

		PUT_TO_JSON(physicsParams, c_GD);
		PUT_TO_JSON(physicsParams, tau_GD);

		PUT_TO_JSON(physicsParams, backtracingLineSearchMaxIters);
		PUT_TO_JSON(physicsParams, backtracingLineSearchAlpha);
		PUT_TO_JSON(physicsParams, backtracingLineSearchTau);
		PUT_TO_JSON(physicsParams, backtracingLineSearchC);

		PUT_TO_JSON(physicsParams, useDouble3x3);
		PUT_TO_JSON(physicsParams, frictionStartIter);

		PUT_TO_JSON(physicsParams, useGPU);
		PUT_TO_JSON(physicsParams, useGDSolver);
		PUT_TO_JSON(physicsParams, useNewton);
		PUT_TO_JSON(physicsParams, useLineSearch);
		PUT_TO_JSON(physicsParams, lineSearchGapIter);
		PUT_TO_JSON(physicsParams, useAccelerator);
		PUT_TO_JSON(physicsParams, acceleratorRho);
		PUT_TO_JSON(physicsParams, GDSolverUseBlockJacobi);
		PUT_TO_JSON(physicsParams, NewtonUseCG);

		// collision detecion
		PUT_TO_JSON(physicsParams, handleCollision);
		PUT_TO_JSON(physicsParams, collisionStiffness);
		PUT_TO_JSON(physicsParams, contactDetectionIters);
		PUT_TO_JSON(physicsParams, collisionEnergyType);
		PUT_TO_JSON(physicsParams, collisionAirDistance);
		PUT_TO_JSON(physicsParams, boundaryCollisionStiffness);
		PUT_TO_JSON(physicsParams, collisionSolutionType);
		PUT_TO_JSON(physicsParams, intermediateCollisionIterations);
		PUT_TO_JSON(physicsParams, boundaryFrictionEpsV);
		PUT_TO_JSON(physicsParams, activeCollisionListPreAllocationRatio);
		PUT_TO_JSON(physicsParams, evaluateConvergence);
		PUT_TO_JSON(physicsParams, evaluationSteps);
		PUT_TO_JSON(physicsParams, saveConvergenceEvaluationResults);

		// for debug
		PUT_TO_JSON(physicsParams, saveIntermediateResults);
		PUT_TO_JSON(physicsParams, debugOutputInterval);
		PUT_TO_JSON(physicsParams, collisionOffHeight);

		PUT_TO_JSON(physicsParams, numThreadsVBDSolve);
		return true;
	}
}
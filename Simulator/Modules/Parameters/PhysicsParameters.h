#pragma once

#include "MeshFrame/Utility/Parser.h"
#include "../Types/Types.h"

namespace GAIA {
	using MF::parseJsonParameters;

	class BasePhysicsParams : public MF::BaseJsonConfig {
	public:
		typedef std::shared_ptr<BasePhysicsParams> BaseSharedPtr;
		typedef BasePhysicsParams* BasePtr;

		typedef std::shared_ptr<BasePhysicsParams> SharedPtr;
		typedef BasePhysicsParams* Ptr;

		BasePhysicsParams() {
			bowlCenter << 0, 0, 0;
			gravity << 0, -10, 0;

			worldBounds << -1e7f, 1e7f,
				-1e7f, 1e7f,
				-1e7f, 1e7f;
		}

		// basics
		Vec3 gravity;

		// step & iterations
		int numFrames = 5000;
		FloatingType timeStep = 0.01666666666666666666;
		int numSubsteps = 20;
		// this is derived from timeStep and numSubsteps, not input parameters
		FloatingType dt = -1.f;
		FloatingType dtSqrReciprocal = -1.f;
		int iterations = 1;

		// Boundary setting
		bool checkAndUpdateWorldBounds = false;
		bool useBowlGround = false;
		bool bowlCap = false;
		bool usePlaneGround = true;
		FloatingType bowlRadius = 1;
		FloatingType bowlOuterRadius = 1.05;
		Vec3 bowlCenter;
		Eigen::Matrix<FloatingType, 3, 2> worldBounds;

		FloatingType boundaryFrictionDynamic = 0.1;
		FloatingType boundaryFrictionStatic = 0.1;

		// debug & outputs
		bool debug;
		int debugVerboseLvl = 1;
		bool showSubstepProgress = false;
		bool showTimeConsumption = false;
		bool outputIntermediateState = false;
		int minParallizationBatchSize = 100;
		std::string debugOutFolder;

		// evaluation & statistics
		bool evaluateConvergence = true;
		bool doStatistics = false;
		bool outputStatistics = false;

		// Parallelization
		bool perMeshParallelization = false;
		bool perTetParallelization = true;

		// Damping
		bool stepInvariantVelDamping = false;
		FloatingType constantVelDamping = -1;

		// Collision Detection
		int collisionDetectionSubSteps = 1;
		bool doCollDetectionOnlyForFirstIteration = true;

		// Collision Solving
		bool smoothSurfaceNormal = true;

		// BVH Updates
		int dcdTetMeshSceneBVHRebuildSteps = 16;
		int dcdSurfaceSceneBVHRebuildSteps = 3;
		int ccdBVHRebuildSteps = 7;

		// input/ouptut
		// master switch
		bool saveOutputs = true;
		bool saveAllModelsTogether = true;
		// output the vertex coordinates (pos only)
		std::string outputExt = "ply";
		int binaryModeVisualizationSteps = 100;

		bool outputT = false;
		bool outputVTK = false;
		int outputStatisticsStep = 10;
		bool outputRecoveryState = true;
		bool saveSimulationParameters = true;
		int outputRecoveryStateStep = 10;
		bool outputRecoveryStateBinary = true;

		// rendering
		std::string shaderFolderPath;

		// newton
		bool psd = true;

		bool fromJson(nlohmann::json& physicsParam) {

			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(physicsParam, gravity);
			PARSE_VEC3_INDEXED_BY_ROUND_BRACKET(physicsParam, bowlCenter);

			std::array<std::array<double, 3>, 2> worldBoundsArr;
			if (parseJsonParameters(physicsParam, "worldBounds", worldBoundsArr))
			{

				worldBounds(0, 0) = worldBoundsArr[0][0];
				worldBounds(1, 0) = worldBoundsArr[0][1];
				worldBounds(2, 0) = worldBoundsArr[0][2];

				worldBounds(0, 1) = worldBoundsArr[1][0];
				worldBounds(1, 1) = worldBoundsArr[1][1];
				worldBounds(2, 1) = worldBoundsArr[1][2];
			}

			EXTRACT_FROM_JSON(physicsParam, numFrames);
			EXTRACT_FROM_JSON(physicsParam, timeStep);
			EXTRACT_FROM_JSON(physicsParam, numSubsteps);
			EXTRACT_FROM_JSON(physicsParam, iterations);
			EXTRACT_FROM_JSON(physicsParam, collisionDetectionSubSteps);
			EXTRACT_FROM_JSON(physicsParam, useBowlGround);
			EXTRACT_FROM_JSON(physicsParam, bowlRadius);
			EXTRACT_FROM_JSON(physicsParam, usePlaneGround);
			EXTRACT_FROM_JSON(physicsParam, boundaryFrictionDynamic);
			EXTRACT_FROM_JSON(physicsParam, boundaryFrictionStatic);


			EXTRACT_FROM_JSON(physicsParam, showSubstepProgress);
			EXTRACT_FROM_JSON(physicsParam, debug);
			EXTRACT_FROM_JSON(physicsParam, debugVerboseLvl);
			EXTRACT_FROM_JSON(physicsParam, showTimeConsumption);
			EXTRACT_FROM_JSON(physicsParam, doCollDetectionOnlyForFirstIteration);
			EXTRACT_FROM_JSON(physicsParam, outputIntermediateState);
			EXTRACT_FROM_JSON(physicsParam, doStatistics);
			EXTRACT_FROM_JSON(physicsParam, debugOutFolder);

			EXTRACT_FROM_JSON(physicsParam, perMeshParallelization);
			EXTRACT_FROM_JSON(physicsParam, perTetParallelization);
			EXTRACT_FROM_JSON(physicsParam, bowlCap);
			EXTRACT_FROM_JSON(physicsParam, bowlOuterRadius);
			EXTRACT_FROM_JSON(physicsParam, smoothSurfaceNormal);
			EXTRACT_FROM_JSON(physicsParam, stepInvariantVelDamping);
			EXTRACT_FROM_JSON(physicsParam, constantVelDamping);

			EXTRACT_FROM_JSON(physicsParam, evaluateConvergence);

			EXTRACT_FROM_JSON(physicsParam, dcdTetMeshSceneBVHRebuildSteps);
			EXTRACT_FROM_JSON(physicsParam, dcdSurfaceSceneBVHRebuildSteps);
			EXTRACT_FROM_JSON(physicsParam, ccdBVHRebuildSteps);

			EXTRACT_FROM_JSON(physicsParam, checkAndUpdateWorldBounds);

			EXTRACT_FROM_JSON(physicsParam, saveOutputs);
			EXTRACT_FROM_JSON(physicsParam, saveAllModelsTogether);
			EXTRACT_FROM_JSON(physicsParam, outputExt);
			EXTRACT_FROM_JSON(physicsParam, shaderFolderPath);
			EXTRACT_FROM_JSON(physicsParam, outputRecoveryState);
			EXTRACT_FROM_JSON(physicsParam, outputRecoveryStateStep);
			EXTRACT_FROM_JSON(physicsParam, outputRecoveryStateBinary);
			EXTRACT_FROM_JSON(physicsParam, outputVTK);
			EXTRACT_FROM_JSON(physicsParam, outputT);
			EXTRACT_FROM_JSON(physicsParam, saveSimulationParameters);
			EXTRACT_FROM_JSON(physicsParam, outputStatistics);
			EXTRACT_FROM_JSON(physicsParam, outputStatisticsStep);
			EXTRACT_FROM_JSON(physicsParam, binaryModeVisualizationSteps);

			EXTRACT_FROM_JSON(physicsParam, psd);
			return true;
		};

		bool toJson(nlohmann::json& physicsParam) {
			PUT_TO_JSON_VEC3_INDEXED_BY_ROUND_BRACKET(physicsParam, gravity);
			PUT_TO_JSON_VEC3_INDEXED_BY_ROUND_BRACKET(physicsParam, bowlCenter);

			std::array<std::array<double, 3>, 2> worldBoundsArr;
			worldBoundsArr[0][0] = worldBounds(0, 0);
			worldBoundsArr[0][1] = worldBounds(1, 0);
			worldBoundsArr[0][2] = worldBounds(2, 0);
			worldBoundsArr[1][0] = worldBounds(0, 1);
			worldBoundsArr[1][1] = worldBounds(1, 1);
			worldBoundsArr[1][2] = worldBounds(2, 1);
			physicsParam["worldBounds"] = worldBoundsArr;

			PUT_TO_JSON(physicsParam, numFrames);
			PUT_TO_JSON(physicsParam, timeStep);
			PUT_TO_JSON(physicsParam, numSubsteps);
			PUT_TO_JSON(physicsParam, iterations);
			PUT_TO_JSON(physicsParam, collisionDetectionSubSteps);
			PUT_TO_JSON(physicsParam, useBowlGround);
			PUT_TO_JSON(physicsParam, bowlRadius);
			PUT_TO_JSON(physicsParam, usePlaneGround);
			PUT_TO_JSON(physicsParam, boundaryFrictionDynamic);
			PUT_TO_JSON(physicsParam, boundaryFrictionStatic);

			PUT_TO_JSON(physicsParam, showSubstepProgress);
			PUT_TO_JSON(physicsParam, debug);
			PUT_TO_JSON(physicsParam, debugVerboseLvl);
			PUT_TO_JSON(physicsParam, showTimeConsumption);
			PUT_TO_JSON(physicsParam, doCollDetectionOnlyForFirstIteration);
			PUT_TO_JSON(physicsParam, outputIntermediateState);
			PUT_TO_JSON(physicsParam, doStatistics);
			PUT_TO_JSON(physicsParam, debugOutFolder);

			PUT_TO_JSON(physicsParam, perMeshParallelization);
			PUT_TO_JSON(physicsParam, perTetParallelization);
			PUT_TO_JSON(physicsParam, bowlCap);
			PUT_TO_JSON(physicsParam, bowlOuterRadius);
			PUT_TO_JSON(physicsParam, smoothSurfaceNormal);
			PUT_TO_JSON(physicsParam, stepInvariantVelDamping);
			PUT_TO_JSON(physicsParam, constantVelDamping);

			PUT_TO_JSON(physicsParam, dcdTetMeshSceneBVHRebuildSteps);
			PUT_TO_JSON(physicsParam, dcdSurfaceSceneBVHRebuildSteps);
			PUT_TO_JSON(physicsParam, ccdBVHRebuildSteps);

			PUT_TO_JSON(physicsParam, evaluateConvergence);
			PUT_TO_JSON(physicsParam, checkAndUpdateWorldBounds);

			PUT_TO_JSON(physicsParam, saveOutputs);
			PUT_TO_JSON(physicsParam, saveAllModelsTogether);
			PUT_TO_JSON(physicsParam, outputExt);
			PUT_TO_JSON(physicsParam, shaderFolderPath);
			PUT_TO_JSON(physicsParam, outputRecoveryState);
			PUT_TO_JSON(physicsParam, outputRecoveryStateStep);
			PUT_TO_JSON(physicsParam, outputRecoveryStateBinary);
			PUT_TO_JSON(physicsParam, outputVTK);
			PUT_TO_JSON(physicsParam, outputT);
			PUT_TO_JSON(physicsParam, saveSimulationParameters);
			PUT_TO_JSON(physicsParam, outputStatistics);
			PUT_TO_JSON(physicsParam, outputStatisticsStep);
			PUT_TO_JSON(physicsParam, binaryModeVisualizationSteps);

			PUT_TO_JSON(physicsParam, psd);

			return true;
		};

	};
}
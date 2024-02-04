#pragma once

#include "../Parameters/PhysicsParameters.h"
#include "PBDTetMeshFEM.h"
namespace GAIA {

	struct PBDPhysicParameters : public BasePhysicsParams {
        /* 
        * material solving
        */

        // controls how many times the inversion solver will attempt to resolve all the verts
        // times of attempts = maxInversionSolveMultiplier * num initially inverted tets
        
        // inversion solving
        bool solveInvertedTets = true;
        FloatingType inversionSolveConstraintMultiplier = 0.02;
        // - for cpu impelmentation
        FloatingType maxInversionSolveMultiplier = 1.5;
        FloatingType  maxSurfaceInversionSolveMultiplier = 2;
        // - for gpu impelmentation
        int inversionSolveIterations = 1;

        int minParallizationBatchSize = 100;
        FloatingType collisionConstraintMultiplier = 1.01;

        // cuda settings
        int numThreadsMaterialSolve = 32;
        int numThreadsInversionSolve = 128;
        int numThreadsBoundarySolve = 128;

        virtual bool fromJson(nlohmann::json& physicsJsonParams) {
            BasePhysicsParams::fromJson(physicsJsonParams);
            EXTRACT_FROM_JSON(physicsJsonParams, solveInvertedTets);
            EXTRACT_FROM_JSON(physicsJsonParams, maxInversionSolveMultiplier);
            EXTRACT_FROM_JSON(physicsJsonParams, maxSurfaceInversionSolveMultiplier);
            EXTRACT_FROM_JSON(physicsJsonParams, minParallizationBatchSize);

            EXTRACT_FROM_JSON(physicsJsonParams, numThreadsMaterialSolve);
            EXTRACT_FROM_JSON(physicsJsonParams, numThreadsBoundarySolve);

            EXTRACT_FROM_JSON(physicsJsonParams, collisionConstraintMultiplier);
            EXTRACT_FROM_JSON(physicsJsonParams, inversionSolveConstraintMultiplier);
            EXTRACT_FROM_JSON(physicsJsonParams, inversionSolveIterations);

            return true;
        };
        virtual bool toJson(nlohmann::json& physicsJsonParams) {
            BasePhysicsParams::toJson(physicsJsonParams);

            PUT_TO_JSON(physicsJsonParams, solveInvertedTets);
            PUT_TO_JSON(physicsJsonParams, maxInversionSolveMultiplier);
            PUT_TO_JSON(physicsJsonParams, maxSurfaceInversionSolveMultiplier);
            PUT_TO_JSON(physicsJsonParams, minParallizationBatchSize);

            PUT_TO_JSON(physicsJsonParams, numThreadsMaterialSolve);
            PUT_TO_JSON(physicsJsonParams, numThreadsBoundarySolve);

            PUT_TO_JSON(physicsJsonParams, collisionConstraintMultiplier);
            PUT_TO_JSON(physicsJsonParams, inversionSolveConstraintMultiplier);
            PUT_TO_JSON(physicsJsonParams, inversionSolveIterations);
            return true;
        };
	};

}
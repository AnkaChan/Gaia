#pragma once

#include "DiscreteCollisionDetector.h"

namespace GAIA {
	struct TetMeshFEM;

	struct ContinuousCollisionDetector {

		ContinuousCollisionDetector(const CollisionDetectionParamters& in_params);
		void initialize(std::vector<std::shared_ptr<TetMeshFEM>> tMeshes);

		bool vertexContinuousCollisionDetection(int32_t vId, int32_t tMeshId, CollisionDetectionResult* pResult);

		void updateBVH(RTCBuildQuality quality);

		std::vector<std::shared_ptr<TetMeshFEM>> tMeshPtrs;
		RTCScene surfaceTriangleTrajectoryScene;
		//RTCScene edgeEdgeTrajectoryScene;
		RTCDevice device;
		const CollisionDetectionParamters& params;
		size_t numFaces;

	};

}
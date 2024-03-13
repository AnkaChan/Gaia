#pragma once
#include "../CollisionDetector/CollisionDetertionParameters.h"
namespace GAIA {

	struct VBDCollisionInfo
	{
		Mat3 frictionHessian;
		Vec3 frictionForce;
		Vec3 diff;
	};

	struct VBDCollisionDetectionResult : public CollisionDetectionResult
	{
		void clear() {
			CollisionDetectionResult::clear();
			collisionForceAndHessian.clear();

		}
		// for collision & friction force and hessian
		CPArray<VBDCollisionInfo, PREALLOCATED_NUM_COLLISIONS> collisionForceAndHessian;

	};

}
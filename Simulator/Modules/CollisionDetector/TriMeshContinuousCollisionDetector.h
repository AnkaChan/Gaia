#pragma once

#include "DiscreteCollisionDetector.h"

#define RECORD_COLLIDING_POINT

namespace GAIA {
	struct TriMeshFEM;

	struct CCDGeometry {
		TriMeshFEM* pMesh = nullptr;
		TVerticesMat* pPrevPos = nullptr;
	};

	struct VFCollision {
		Vec3 barycentrics;
#ifdef RECORD_COLLIDING_POINT
		Vec3 c;
#endif // RECORD_COLLIDING_POINT

		FloatingType t;

		int vertexId;
		int vertexMeshId;

		int faceId;
		int faceMeshId;
	};

	struct EECollision {
#ifdef RECORD_COLLIDING_POINT
		Vec3 c;
#endif // RECORD_COLLIDING_POINT
		FloatingType miu1; // colliding point on edge 1 p = e1.p1 + miu1 * (e1.p2 - e1.p1)
		FloatingType miu2;

		FloatingType t;
	};

	template<typename T>
	struct CollisionResults {
		void initialize(size_t initialMaxNumVFCollisions)
		{
			maxNumTriTriIntersections = initialMaxNumVFCollisions;
			collisions.resize(maxNumTriTriIntersections);
		}

		void clear() {
			numCollisions = 0;
			if (overflow) {
				maxNumTriTriIntersections = size_t(1.5 * collisions.size());
				collisions.resize(maxNumTriTriIntersections);
				overflow = false;
			}
		}

		size_t maxNumTriTriIntersections = -1;
		std::atomic<int> numCollisions = 0;
		bool overflow = false;
		std::vector<T> collisions;

	};



	typedef CollisionResults<VFCollision> VFCollisionResults;
	typedef CollisionResults<EECollision> EECollisionResults;

	struct TriMeshContinuousCollisionDetector {
		typedef std::shared_ptr<TriMeshContinuousCollisionDetector> SharedPtr;
		typedef TriMeshContinuousCollisionDetector* Ptr;

		TriMeshContinuousCollisionDetector(const CollisionDetectionParamters& in_params);
		void initialize(std::vector<std::shared_ptr<TriMeshFEM>> tMeshes, std::vector<TVerticesMat*> prevPoses);

		bool continuousCollisionDetection();

		void updateBVH(RTCBuildQuality quality);

		RTCScene triangleTrajectoriesScene;
		RTCScene vertexTrajectoriesScene;
		RTCScene edgeTrajectoriesScene;
		RTCDevice device;
		const CollisionDetectionParamters& params;
		size_t numFaces;

		std::vector<CCDGeometry> ccdGeometries;
		VFCollisionResults vfCollisionResults;
		EECollisionResults eeCollisionResults;

		float preAllocationRatio=0.25f;
	};

}
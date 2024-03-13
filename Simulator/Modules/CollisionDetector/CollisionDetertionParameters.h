#pragma once

#include "MeshFrame/Utility/Parser.h"
#include <MeshFrame/Memory/Array.h>
#define PREALLOCATED_NUM_COLLISIONS 1

#include "../Types/Types.h"

namespace GAIA {

	struct CollisionStatistics : public MF::BaseJsonConfig
	{
		// number of total sub-timesteps
		std::vector<int> numOfCollisionsDCDs;
		std::vector<int> numOfCollisionsCCDs;
		// number of total sub-timesteps x number of DCD result
		std::vector<std::vector<int>> numOfBVHQuerysEachStep;
		std::vector<std::vector<int>> numOfTetTraversal;
		std::vector<std::vector<int>> numTetsTraversed;

		virtual bool fromJson(nlohmann::json& collisionParam) {
			EXTRACT_FROM_JSON(collisionParam, numOfCollisionsDCDs);
			EXTRACT_FROM_JSON(collisionParam, numOfCollisionsCCDs);
			EXTRACT_FROM_JSON(collisionParam, numOfBVHQuerysEachStep);
			EXTRACT_FROM_JSON(collisionParam, numOfTetTraversal);
			EXTRACT_FROM_JSON(collisionParam, numTetsTraversed);

			return true;
		}

		virtual bool toJson(nlohmann::json& collisionParam) {
			PUT_TO_JSON(collisionParam, numOfCollisionsDCDs);
			PUT_TO_JSON(collisionParam, numOfCollisionsCCDs);
			PUT_TO_JSON(collisionParam, numOfBVHQuerysEachStep);
			PUT_TO_JSON(collisionParam, numOfTetTraversal);
			PUT_TO_JSON(collisionParam, numTetsTraversed);

			return true;

		}
	};

	struct VolumetricCollisionParameters : public MF::BaseJsonConfig {
		int numPreallocatedTritriIntersections = 262144 / 4;
		float mergeTolerance = 1e-4f;
		float rayPertubation = 2e-2f;

		virtual bool fromJson(nlohmann::json& physicsJsonParams) {
			EXTRACT_FROM_JSON(physicsJsonParams, numPreallocatedTritriIntersections);
			EXTRACT_FROM_JSON(physicsJsonParams, mergeTolerance);

			return true;
		};
		virtual bool toJson(nlohmann::json& physicsJsonParams) {
			PUT_TO_JSON(physicsJsonParams, numPreallocatedTritriIntersections);
			PUT_TO_JSON(physicsJsonParams, mergeTolerance);
			return true;
		}
	};

	struct CollisionDetectionParamters : public MF::BaseJsonConfig
	{
		// collision detection parameters
		bool allowCCD = true;
		bool allowDCD = true;

		// DCD parameters
		bool checkFeasibleRegion = true;
		bool checkTetTraverse = true;
		bool handleSelfCollision = true;
		bool stopTraversingAfterPassingQueryPoint = true;
		bool tetrahedralTraverseForNonSelfIntersection = true;
		bool useStaticTraverse = true;
		bool restPoseCloestPoint = false;
		bool loopLessTraverse = false;
		int maxNumberOfCollisions = 1;

		// CCD parameters
		bool doEdgeEdgeCCD = false;

		bool shiftQueryPointToCenter = true;
		float centerShiftLevel = 0.01f;
		int maxNumberOfBVHQuery = 500;

		// collision info computation
		bool computeContactNormal = false;

		//int numberOfBVHQuery = 0;
		//int numberOfTetTraversal = 0;
		//int numberOfTetsTraversed = 0;

		// tetrahedral traverse parameters
		float feasibleRegionEpsilon = 1e-2f;
		float rayTriIntersectionEPSILON = 1e-10f;
		float maxSearchDistanceMultiplier = 1.8f;

		// volumetric related collision
		bool allowVolumetricCollision = false;
		VolumetricCollisionParameters volCollisionParams;

		virtual bool fromJson(nlohmann::json& collisionParam) {
			EXTRACT_FROM_JSON(collisionParam, allowCCD);
			EXTRACT_FROM_JSON(collisionParam, allowDCD);
			EXTRACT_FROM_JSON(collisionParam, checkTetTraverse);
			EXTRACT_FROM_JSON(collisionParam, checkFeasibleRegion);
			EXTRACT_FROM_JSON(collisionParam, tetrahedralTraverseForNonSelfIntersection);
			EXTRACT_FROM_JSON(collisionParam, handleSelfCollision);
			EXTRACT_FROM_JSON(collisionParam, stopTraversingAfterPassingQueryPoint);
			EXTRACT_FROM_JSON(collisionParam, maxNumberOfCollisions);

			EXTRACT_FROM_JSON(collisionParam, computeContactNormal);

			EXTRACT_FROM_JSON(collisionParam, shiftQueryPointToCenter);
			EXTRACT_FROM_JSON(collisionParam, centerShiftLevel);

			EXTRACT_FROM_JSON(collisionParam, feasibleRegionEpsilon);
			EXTRACT_FROM_JSON(collisionParam, rayTriIntersectionEPSILON);
			EXTRACT_FROM_JSON(collisionParam, maxSearchDistanceMultiplier);
			EXTRACT_FROM_JSON(collisionParam, useStaticTraverse);

			EXTRACT_FROM_JSON(collisionParam, restPoseCloestPoint);
			EXTRACT_FROM_JSON(collisionParam, loopLessTraverse);

			EXTRACT_FROM_JSON(collisionParam, allowVolumetricCollision);


			return true;
		}

		virtual bool toJson(nlohmann::json& collisionParam) {
			PUT_TO_JSON(collisionParam, allowCCD);
			PUT_TO_JSON(collisionParam, allowDCD);
			PUT_TO_JSON(collisionParam, checkTetTraverse);
			PUT_TO_JSON(collisionParam, checkFeasibleRegion);
			PUT_TO_JSON(collisionParam, tetrahedralTraverseForNonSelfIntersection);
			PUT_TO_JSON(collisionParam, handleSelfCollision);
			PUT_TO_JSON(collisionParam, stopTraversingAfterPassingQueryPoint);
			PUT_TO_JSON(collisionParam, maxNumberOfCollisions);

			PUT_TO_JSON(collisionParam, computeContactNormal);

			PUT_TO_JSON(collisionParam, shiftQueryPointToCenter);
			PUT_TO_JSON(collisionParam, centerShiftLevel);

			PUT_TO_JSON(collisionParam, feasibleRegionEpsilon);
			PUT_TO_JSON(collisionParam, rayTriIntersectionEPSILON);
			PUT_TO_JSON(collisionParam, maxSearchDistanceMultiplier);
			PUT_TO_JSON(collisionParam, useStaticTraverse);

			PUT_TO_JSON(collisionParam, restPoseCloestPoint);
			PUT_TO_JSON(collisionParam, loopLessTraverse);

			PUT_TO_JSON(collisionParam, allowVolumetricCollision);

			return true;

		}
	};

	enum class ClosestPointOnTriangleType
	{
		AtA,
		AtB,
		AtC,
		AtAB,
		AtBC,
		AtAC,
		AtInterior,
		NotFound
	};

	enum class ClosestPointOnPrimitiveType
	{
		AtVertex,
		AtEdge,
		AtFace,
		AtInterior,
		NotFound
	};

	ClosestPointOnPrimitiveType getClosestPointOnPrimitiveType(const ClosestPointOnTriangleType& type);

	struct DiscreteCollisionDetector;

	struct CollidingPointInfo {
		// for DCD only
		Vec3 closestSurfacePt;
		Vec3 closestSurfacePtBarycentrics;
		// for DCD + CCD
		Vec3 closestPointNormal;

		IdType intersectedElement;
		IdType intersectedMeshId;

		IdType closestSurfaceFaceId;
		ClosestPointOnTriangleType closestPointType;

		bool shortestPathFound;
	};

	struct CollisionDetectionResult
	{
		CollisionDetectionResult()
		{}

		int numIntersections() { return collidingPts.size(); }
		void clear() {
			collidingPts.clear();

			//numberOfBVHQuery = 0;
			//numberOfTetTraversal = 0;
			//numberOfTetsTraversed = 0;
		}

		CPArray<CollidingPointInfo, PREALLOCATED_NUM_COLLISIONS> collidingPts;


		// set to non-negative when doing vertex collision detection
		int idVQuery = -1;
		// set to non-nullptr when doing tet centroid collision detection
		int idTetQuery = -1;

		int idTMQuery = -1;

		//std::map<unsigned int, PathFinder::TM::Ptr>* pTetmeshGeoIdToPointerMap;
		//std::map<PathFinder::TM::Ptr, unsigned int>* pTetmeshPtrToTetMeshIndexMap;
		//std::vector<std::vector<PathFinder::TM::TPtr>>* pTetTraversed = nullptr;
		//std::vector<PathFinder::TM::Ptr>* pTetMeshPtrs;
		bool handleSelfIntersection = true;
		bool fromCCD = false;

		float penetrationDepth = -1.f;


		//int numberOfBVHQuery = 0;
		//int numberOfTetTraversal = 0;
		//int numberOfTetsTraversed = 0;

		// either CCD or DCD
		void* pDetector = nullptr;
	};

	inline ClosestPointOnPrimitiveType GAIA::getClosestPointOnPrimitiveType(const ClosestPointOnTriangleType& type)
	{
		switch (type)
		{
		case ClosestPointOnTriangleType::AtA:
		case ClosestPointOnTriangleType::AtB:
		case ClosestPointOnTriangleType::AtC:
			return ClosestPointOnPrimitiveType::AtVertex;
		case ClosestPointOnTriangleType::AtAB:
		case ClosestPointOnTriangleType::AtBC:
		case ClosestPointOnTriangleType::AtAC:
			return ClosestPointOnPrimitiveType::AtEdge;
		case ClosestPointOnTriangleType::AtInterior:
			return ClosestPointOnPrimitiveType::AtFace;
		case ClosestPointOnTriangleType::NotFound:
			return ClosestPointOnPrimitiveType::NotFound;
		default:
			return ClosestPointOnPrimitiveType::NotFound;
		}
	}
}
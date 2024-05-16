#pragma once
#include "../TriMesh/TriMesh.h"

namespace GAIA {
	struct ColliderTrimeshBaseParams : public TriMeshParams
	{
		typedef std::shared_ptr<ColliderTrimeshBaseParams> BaseSharedPtr;
		typedef ColliderTrimeshBaseParams* BasePtr;
		typedef std::shared_ptr<ColliderTrimeshBaseParams> SharedPtr;
		typedef ColliderTrimeshBaseParams* Ptr;

		std::string colliderType;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			TriMeshParams::fromJson(objectParam);
			EXTRACT_FROM_JSON(objectParam, colliderType);
			return true;
		}

		inline bool toJson(nlohmann::json& objectParam)
		{
			TriMeshParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, colliderType);
			return true;
		}
	};

	struct ColliderTriMeshBase : public TriMeshFEM
	{
		typedef std::shared_ptr<ColliderTriMeshBase> SharedPtr;
		typedef ColliderTriMeshBase* Ptr;

		virtual void update(IdType frameId, IdType substepId, IdType iter, size_t numSubsteps, size_t numIters) = 0;
		virtual void initialize(ColliderTrimeshBaseParams::SharedPtr inObjectParams)
		{
			pParams = inObjectParams;
		};

		ColliderTrimeshBaseParams::SharedPtr pParams;
	};
}
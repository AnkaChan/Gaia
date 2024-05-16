#pragma once
#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"

namespace GAIA {
	struct ColliderTetMeshBaseParams : public ObjectParams
	{
		typedef std::shared_ptr<ColliderTetMeshBaseParams> BaseSharedPtr;
		typedef ColliderTetMeshBaseParams* BasePtr;
		typedef std::shared_ptr<ColliderTetMeshBaseParams> SharedPtr;
		typedef ColliderTetMeshBaseParams* Ptr;

		std::string colliderType;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			ObjectParams::fromJson(objectParam);
			EXTRACT_FROM_JSON(objectParam, colliderType);
			return true;
		}

		inline bool toJson(nlohmann::json& objectParam)
		{
			ObjectParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, colliderType);
			return true;
		}
	};

	struct ColliderTetMeshBase : public TetMeshFEM
	{
		typedef std::shared_ptr<ColliderTetMeshBase> SharedPtr;
		typedef ColliderTetMeshBase* Ptr;

		virtual void update(IdType frameId, IdType substepId, IdType iter, size_t numSubsteps, size_t numIters) = 0;
		virtual void initialize(ColliderTetMeshBaseParams::SharedPtr inObjectParams)
		{
			pParams = inObjectParams;
		};

		virtual void initilizeSurfaceTriMesh();
		virtual void updateSurfaceTriMesh();

		TriMeshFEM::SharedPtr pSurfaceTriMesh;

		ColliderTetMeshBaseParams::SharedPtr pParams;
	};
}
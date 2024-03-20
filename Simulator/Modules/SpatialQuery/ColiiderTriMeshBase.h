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
	};

	struct ColliderTrimeshBase : public TriMeshFEM
	{
		typedef std::shared_ptr<ColliderTrimeshBase> SharedPtr;
		typedef ColliderTrimeshBase* Ptr;

		virtual void update(IdType frameId, IdType substepId, IdType iter) = 0;
		virtual void initialize(ColliderTrimeshBaseParams::SharedPtr inObjectParams)
		{
			pParams = inObjectParams;
		};

		ColliderTrimeshBaseParams::SharedPtr pParams;
	};
}
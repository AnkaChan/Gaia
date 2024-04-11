#include "ColiiderTriMeshBase.h"

namespace GAIA {
	struct ColliderTrimeshSequenceParams : ColliderTrimeshBaseParams
	{
		typedef std::shared_ptr<ColliderTrimeshSequenceParams> SharedPtr;
		typedef ColliderTrimeshSequenceParams* Ptr;
		std::vector<std::string> meshFiles;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			ColliderTrimeshBaseParams::fromJson(objectParam);
			EXTRACT_FROM_JSON(objectParam, meshFiles);
			return true;
		}

		inline bool toJson(nlohmann::json& objectParam)
		{
			ColliderTrimeshBaseParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, meshFiles);
			return true;
		}

	};

	// the simpliest type of collider, which is a sequence of mesh files
	// currently, only mesh files share the same topology are supported
	struct ColliderTrimeshSequence : public ColliderTrimeshBase
	{
		virtual void update(IdType frameId, IdType substepId, IdType iter) 
		{
			if (frameId != curFrameId && frameId >=0 && frameId < colliderParameters().meshFiles.size())
			{
				loadObj(colliderParameters().meshFiles[frameId]);
				updated = true;
			}
			else
			{
				updated = false;
			}
		};
		virtual void initialize(ColliderTrimeshBaseParams::SharedPtr inObjectParams) 
		{
			ColliderTrimeshBase::initialize(inObjectParams);
			pParams = inObjectParams;
			// no need to call the base class' initialization function, because it's not used for simulation
			if (colliderParameters().meshFiles.size())
			{
				inObjectParams->path = colliderParameters().meshFiles[0];
				TriMeshFEM::initialize(inObjectParams, true);
				curFrameId = 0;
			}
		};

		ColliderTrimeshSequenceParams& colliderParameters() {
			return *(ColliderTrimeshSequenceParams*)pParams.get();
		}

		IdType curFrameId = -1;
	};
}
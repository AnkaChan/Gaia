#include "ColiiderTriMeshBase.h"

namespace GAIA {
	struct ColliderTrimeshSequenceParams : ColliderTrimeshBaseParams
	{
		typedef std::shared_ptr<ColliderTrimeshSequenceParams> SharedPtr;
		typedef ColliderTrimeshSequenceParams* Ptr;
		std::vector<std::string> meshFiles;

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
			}
		};
		virtual void initialize(ColliderTrimeshBaseParams::SharedPtr inObjectParams, bool precomputeToplogy) 
		{
			pParams = inObjectParams;
			// no need to call the base class' initialization function, because it's not used for simulation
			if (colliderParameters().meshFiles.size())
			{
				loadObj(colliderParameters().meshFiles[0]);
				curFrameId = 0;
			}
		};

		ColliderTrimeshSequenceParams& colliderParameters() {
			return *(ColliderTrimeshSequenceParams*)pParams.get();
		}

		IdType curFrameId = -1;
	};
}
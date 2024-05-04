#include "ColiiderTriMeshBase.h"

namespace GAIA {
	struct ColliderTrimeshSequenceParams : ColliderTrimeshBaseParams
	{
		typedef std::shared_ptr<ColliderTrimeshSequenceParams> SharedPtr;
		typedef ColliderTrimeshSequenceParams* Ptr;
		std::vector<std::string> meshFiles;
		bool interpolate = true;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			ColliderTrimeshBaseParams::fromJson(objectParam);
			EXTRACT_FROM_JSON(objectParam, meshFiles);
			EXTRACT_FROM_JSON(objectParam, interpolate);
			return true;
		}

		inline bool toJson(nlohmann::json& objectParam)
		{
			ColliderTrimeshBaseParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, meshFiles);
			PUT_TO_JSON(objectParam, interpolate);
			return true;
		}

	};

	// the simpliest type of collider, which is a sequence of mesh files
	// currently, only mesh files share the same topology are supported
	struct ColliderTrimeshSequence : public ColliderTrimeshBase
	{
		virtual void update(IdType frameId, IdType substepId, IdType iter, size_t numsubsteps, size_t numIters) 
		{
			if (frameId != curFrameId && frameId >=0 && frameId+1 < colliderParameters().meshFiles.size())
			{
				curFrameId = frameId;
				curFrameMesh.positions() = nextFrameMesh.positions();
				nextFrameMesh.loadObj(colliderParameters().meshFiles[frameId + 1]);
			}

			if (iter == 0)
			{
				if (colliderParameters().interpolate)
				{
					FloatingType t = FloatingType(numsubsteps - substepId) / numsubsteps;
					positions() = curFrameMesh.positions() * (1 - t) + nextFrameMesh.positions() * (t);
				}
				else
				{
					positions() = curFrameMesh.positions();
				}
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
				curFrameMesh.loadObj(colliderParameters().meshFiles[0]);

				curFrameId = 0;
				if (colliderParameters().meshFiles.size()>=2)
				{
					nextFrameMesh.loadObj(colliderParameters().meshFiles[1]);
				}
				else
				{
					nextFrameMesh.loadObj(colliderParameters().meshFiles[0]);
				}
			}
		};

		ColliderTrimeshSequenceParams& colliderParameters() {
			return *(ColliderTrimeshSequenceParams*)pParams.get();
		}

		IdType curFrameId = -1;
		TriMeshFEM curFrameMesh;
		TriMeshFEM nextFrameMesh;
	};
}
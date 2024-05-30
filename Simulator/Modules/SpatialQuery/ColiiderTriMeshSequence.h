#include "ColiiderTriMeshBase.h"

namespace GAIA {
	struct ColliderTrimeshSequenceParams : ColliderTriMeshBaseParams
	{
		typedef std::shared_ptr<ColliderTrimeshSequenceParams> SharedPtr;
		typedef ColliderTrimeshSequenceParams* Ptr;
		std::vector<std::string> meshFiles{};
		std::vector<int> keyFrames{};
		bool interpolateSubstep = true;
		int interpolateIter = 1;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			ColliderTriMeshBaseParams::fromJson(objectParam);
			EXTRACT_FROM_JSON(objectParam, meshFiles);
			EXTRACT_FROM_JSON(objectParam, interpolateSubstep);
			EXTRACT_FROM_JSON(objectParam, interpolateIter);
			EXTRACT_FROM_JSON(objectParam, keyFrames);
			assert(interpolateIter >= 1);
			if (interpolateIter > 1 && !interpolateSubstep) {
				interpolateSubstep = true;
				std::cout << "Warning: interpolateIter is set to be larger than 0, but interpolateSubstep is set to be false, set interpolateSubstep to be true" << std::endl;
			}
			if (keyFrames.size() == 0)
			{
				keyFrames.resize(meshFiles.size());
				std::iota(keyFrames.begin(), keyFrames.end(), 0);
			}
			if (keyFrames.size() > meshFiles.size()) {
				keyFrames.resize(meshFiles.size());
				std::cout << "Warning: keyFrames size is larger than meshFiles size, resize to meshFiles size" << std::endl;
			}
			for (int i = 1; i < keyFrames.size(); ++i) {
				if (keyFrames[i] <= keyFrames[i - 1]) {
					std::cout << "Error: keyFrames should be in increasing order" << std::endl;
					return false;
				}
			}
			return true;
		}

		inline bool toJson(nlohmann::json& objectParam)
		{
			ColliderTriMeshBaseParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, meshFiles);
			PUT_TO_JSON(objectParam, interpolateSubstep);
			PUT_TO_JSON(objectParam, interpolateIter);
			PUT_TO_JSON(objectParam, keyFrames);
			return true;
		}

	};

	// the simpliest type of collider, which is a sequence of mesh files
	// currently, only mesh files share the same topology are supported
	struct ColliderTrimeshSequence : public ColliderTriMeshBase
	{
		virtual void update(IdType frameId, IdType substepId, IdType iter, size_t numsubsteps, size_t numIters)
		{
			assert(colliderParameters().interpolateIter < numIters);
			// do not update during iterations, only update at the first iteration
			if (frameId < colliderParameters().keyFrames.front() || frameId >= colliderParameters().keyFrames.back())
			{
				updated = false;
				return;
			}
			auto pos = std::upper_bound(colliderParameters().keyFrames.begin(), colliderParameters().keyFrames.end(), frameId);
			const auto prevFrameId = *(pos - 1);
			const auto nextFrameId = *pos;
			const auto& prevPos = meshes[std::distance(colliderParameters().keyFrames.begin(), pos) - 1].positions();
			const auto& nextPos = meshes[std::distance(colliderParameters().keyFrames.begin(), pos)].positions();
			int numFrames = nextFrameId - prevFrameId;
			if (colliderParameters().interpolateSubstep)
			{
				if (iter < colliderParameters().interpolateIter)
				{
					FloatingType t = FloatingType(frameId - prevFrameId) / numFrames;
					t += FloatingType(substepId) / (numsubsteps * numFrames);
					t += FloatingType(iter) / (colliderParameters().interpolateIter * numsubsteps * numFrames);
					positions() = prevPos * (1 - t) + nextPos * t;
					updated = true;
				}
				else
				{
					updated = false;
				}
			}
			else
			{
				if (substepId == 0 && iter == 0) {
					FloatingType t = FloatingType(frameId - prevFrameId) / numFrames;
					positions() = prevPos * (1 - t) + nextPos * t;
					updated = true;
				}
				else {
					updated = false;
				}
			}
		};
		virtual void initialize(ColliderTriMeshBaseParams::SharedPtr inObjectParams)
		{
			ColliderTriMeshBase::initialize(inObjectParams);
			pParams = inObjectParams;
			// no need to call the base class' initialization function, because it's not used for simulation
			if (colliderParameters().keyFrames.size())
			{
				meshes.resize(colliderParameters().keyFrames.size());
				inObjectParams->path = colliderParameters().meshFiles[0];
				TriMeshFEM::initialize(inObjectParams, true);
				for (int i = 0; i < colliderParameters().keyFrames.size(); ++i) {
					inObjectParams->path = colliderParameters().meshFiles[i];
					meshes[i].pObjectParams = inObjectParams;
					meshes[i].loadObj(inObjectParams->path);
					meshes[i].applyRotationScalingTranslation();
				}
			}
		};

		ColliderTrimeshSequenceParams& colliderParameters() {
			return *(ColliderTrimeshSequenceParams*)pParams.get();
		}
		std::vector<TriMeshFEM> meshes{};
	};
}
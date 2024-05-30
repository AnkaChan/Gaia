#pragma once

#include "ColliderTetMeshBase.h"

namespace GAIA
{
	struct ColliderTetMeshSequenceParams : public ColliderTetMeshBaseParams
	{
		typedef std::shared_ptr<ColliderTetMeshSequenceParams> SharedPtr;
		typedef ColliderTetMeshSequenceParams* Ptr;
		std::vector<std::string> meshFiles{};
		std::vector<int> keyFrames{};
		bool interpolateSubstep = true;
		int interpolateIter = 1;

		inline bool fromJson(nlohmann::json& objectParam)
		{
			ColliderTrimeshBaseParams::fromJson(objectParam);
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
			ColliderTrimeshBaseParams::toJson(objectParam);
			PUT_TO_JSON(objectParam, meshFiles);
			PUT_TO_JSON(objectParam, interpolateSubstep);
			PUT_TO_JSON(objectParam, interpolateIter);
			PUT_TO_JSON(objectParam, keyFrames);
			return true;
		}

	};
}
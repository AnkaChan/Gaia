#pragma once

#include <Types/Types.h>

namespace GAIA {
	struct VertexFaceAttachmentConstraintParams : MF::BaseJsonConfig {
		// constraints
		std::vector<IdType> vertexConstraintVertexInfos;        // (meshId, vertex)*N
		std::vector<IdType> vertexConstraintTargetFacesInfos;   // (meshId, FV1, FV2, VF3)*N
		std::vector<std::array<FloatingType, 3>> vertexConstraintBarys;
		std::vector<FloatingType> vertexConstraintStiffness;

		virtual bool fromJson(nlohmann::json& objectParam) {
			EXTRACT_FROM_JSON(objectParam, vertexConstraintVertexInfos);
			EXTRACT_FROM_JSON(objectParam, vertexConstraintTargetFacesInfos);
			EXTRACT_FROM_JSON(objectParam, vertexConstraintBarys);
			EXTRACT_FROM_JSON(objectParam, vertexConstraintStiffness);
			return true;
		};
		virtual bool toJson(nlohmann::json& objectParam) {
			PUT_TO_JSON(objectParam, vertexConstraintVertexInfos);
			PUT_TO_JSON(objectParam, vertexConstraintTargetFacesInfos);
			PUT_TO_JSON(objectParam, vertexConstraintBarys);
			PUT_TO_JSON(objectParam, vertexConstraintStiffness);
			return true;
		}


	};
	
	struct VertexFaceAttachmentConstraint {
		FloatingType constraintStiffness = 1e4;
		
		IdType vertexOrder = -1; // 0~2 face side; 1: vertex Side
		
		IdType vertexId = -1;
		IdType attachedMeshVSdie = -1; 

		IdType attachedMeshFaceSide = -1; // =-1 fixed to a position; >=0 : fixed to a mesh
		Vec3I attachedFaceVIds;
		Vec3 attachedPos;        // in case of 0 fixed: 1: fixed to a mesh
	};




	struct VertexConstraints {
		std::vector<VertexFaceAttachmentConstraint> attachmentConstraints;
	};
}
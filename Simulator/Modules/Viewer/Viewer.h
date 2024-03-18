#pragma once
#include <memory>
#include "../TriMesh/TriMesh.h"
#include "../TetMesh/TetMeshFEM.h"
#include <MeshFrame/Utility/Parser.h>

namespace GAIA {
	struct ViewerParams : public MF::BaseJsonConfig {
		typedef ViewerParams* Ptr;
		typedef std::shared_ptr<ViewerParams> SharedPtr;

		bool enableViewer = false;
		bool enableGround = false;
		float GUIRegularFontSize = 28.f;
		float GUIMonoFontSize = 26.f;

		bool fromJson(nlohmann::json& physicsParam);
		bool toJson(nlohmann::json& physicsParam);


	};

	struct ViewerImpl;

	struct Viewer {
		Viewer(ViewerParams::SharedPtr pViewerParams_in);
		~Viewer();
		void init();
		void update();
		void registerTrimeshes(const std::vector<TriMeshFEM::SharedPtr>& inTrimeshes);
		void registerTetmeshes(const std::vector<TetMeshFEM::SharedPtr>& inTetmeshes);
		
		void frameTick();

		std::shared_ptr<ViewerImpl> pImpl;
		ViewerParams::SharedPtr pViewerParams;
	};
}
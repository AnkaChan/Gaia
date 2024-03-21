#include "Viewer.h"
#include <stdexcept>

void GUINoCompliationError() {
	std::cout << "Error! Trying to call GUI while GUI is not complied!\n"
		<< "You should either remove the \"--gui\" command or turn off the GAIA_NO_GUI compilation option." << std::endl;
	throw std::exception("Trying to call uncompiled GUI components.");
	std::exit(-1);
}


bool GAIA::ViewerParams::fromJson(nlohmann::json& physicsParam)
{
	return false;
}
bool GAIA::ViewerParams::toJson(nlohmann::json& physicsParam)
{
	return false;
}

#ifndef GAIA_NO_GUI
#include <polyscope/polyscope.h>
#include <polyscope/surface_mesh.h>
#include <MeshFrame/Utility/Str.h>
#include <iostream>

namespace polyscope {
	// Forward declare compressed binary font functions
	namespace render {
		unsigned int getCousineRegularCompressedSize();
		const unsigned int* getCousineRegularCompressedData();
		unsigned int getLatoRegularCompressedSize();
		const unsigned int* getLatoRegularCompressedData();
	} // namespace render
}

namespace GAIA {
	struct SetImGUiFontsCallBack{
		SetImGUiFontsCallBack(ViewerParams::SharedPtr pViewerParams_in)
			: pViewerParams(pViewerParams_in)
		{}
		std::tuple<ImFontAtlas*, ImFont*, ImFont*> operator()() 
		{
			ImGuiIO& io = ImGui::GetIO();

			// outputs
			ImFontAtlas* globalFontAtlas;
			ImFont* regularFont;
			ImFont* monoFont;

			{ // add regular font
				ImFontConfig config;
				regularFont = io.Fonts->AddFontFromMemoryCompressedTTF(polyscope::render::getLatoRegularCompressedData(),
					polyscope::render::getLatoRegularCompressedSize(), pViewerParams->GUIRegularFontSize, &config);
			}

			{ // add mono font
				ImFontConfig config;
				monoFont = io.Fonts->AddFontFromMemoryCompressedTTF(polyscope::render::getCousineRegularCompressedData(),
					polyscope::render::getCousineRegularCompressedSize(), pViewerParams->GUIMonoFontSize, &config);
			}

			// io.Fonts->AddFontFromFileTTF("test-font-name.ttf", 16);

			io.Fonts->Build();
			globalFontAtlas = io.Fonts;

			return std::tuple<ImFontAtlas*, ImFont*, ImFont*>{globalFontAtlas, regularFont, monoFont};
		}

		ViewerParams::SharedPtr pViewerParams;
	};


	struct PolyScopeSurfaceMesh
	{
		PolyScopeSurfaceMesh()
			: updated(false)
		{}
		PolyScopeSurfaceMesh(const PolyScopeSurfaceMesh& other) 
			: updated(true)
			, meshType(other.meshType)
			, meshId(other.meshId)
			, pPolyScopeSurfaceMesh(other.pPolyScopeSurfaceMesh)
			, meshName(other.meshName)
		{}
		int meshType = -1; // 0: tetmesh, 1: triMesh
		int meshId = -1;
		polyscope::SurfaceMesh* pPolyScopeSurfaceMesh = nullptr;
		std::atomic<bool> updated ;
		std::string meshName;
	};

	struct ViewerImpl {
		ViewerImpl(ViewerParams::SharedPtr pViewerParams_in) 
			: pViewerParams(pViewerParams_in)
			, setFontcallBack(pViewerParams_in)
		{}
		void update();
		void registerTrimeshes(const std::vector<TriMeshFEM::SharedPtr>& inTrimeshes);
		void registerTetmeshes(const std::vector<TetMeshFEM::SharedPtr>& inTetmeshes);
		void setAllMeshesToUpdated()
		{
			for (size_t iPSMesh = 0; iPSMesh < polyscopeMeshes.size(); iPSMesh++)
			{
				PolyScopeSurfaceMesh& polyScopeSurfaceMesh = polyscopeMeshes[iPSMesh];
				polyScopeSurfaceMesh.updated = true;
			}
		}

		std::vector<TriMeshFEM::SharedPtr> trimeshes;
		std::vector<TetMeshFEM::SharedPtr> tetmeshes;
		std::vector<PolyScopeSurfaceMesh> polyscopeMeshes;

		SetImGUiFontsCallBack setFontcallBack;
		ViewerParams::SharedPtr pViewerParams;
	};

	void ViewerImpl::update()
	{
		for (size_t iPSMesh = 0; iPSMesh < polyscopeMeshes.size(); iPSMesh++)
		{
			PolyScopeSurfaceMesh& polyScopeSurfaceMesh = polyscopeMeshes[iPSMesh];
			if (polyScopeSurfaceMesh.updated)
			{
				polyScopeSurfaceMesh.updated = false;
				GAIA::TriMeshFEM::SharedPtr pTriMesh = nullptr;
				GAIA::TetMeshFEM::SharedPtr pTetMesh = nullptr;

				switch (polyScopeSurfaceMesh.meshType)
				{
				case 0:
					pTetMesh = tetmeshes[polyScopeSurfaceMesh.meshId];
					polyScopeSurfaceMesh.pPolyScopeSurfaceMesh->updateVertexPositions(pTetMesh->positions().transpose());
					break;

				case 1:
					pTriMesh = trimeshes[polyScopeSurfaceMesh.meshId];
					polyScopeSurfaceMesh.pPolyScopeSurfaceMesh->updateVertexPositions(pTriMesh->positions().transpose());
					break;

				default:
					break;
				}
			}
		}
	}

	void ViewerImpl::registerTrimeshes(const std::vector<TriMeshFEM::SharedPtr>& inTrimeshes)
	{
		trimeshes = inTrimeshes;
		for (size_t iMesh = 0; iMesh < trimeshes.size(); iMesh++)
		{
			GAIA::TriMeshFEM::SharedPtr pTriMesh = trimeshes[iMesh];
			std::string meshName = std::to_string(iMesh);
			MF::STR::padTo(meshName, 4, '0');
			meshName = "Trimesh_" + meshName;
			polyscope::SurfaceMesh* pPolyScopeSurfaceMesh =
				polyscope::registerSurfaceMesh(meshName, pTriMesh->positions().transpose(), pTriMesh->facePos.transpose());

			polyscopeMeshes.emplace_back();
			PolyScopeSurfaceMesh& polyScopeSurfaceMesh = polyscopeMeshes.back();
			polyScopeSurfaceMesh.meshName = meshName;
			polyScopeSurfaceMesh.meshType = 1;
			polyScopeSurfaceMesh.meshId = iMesh;
			polyScopeSurfaceMesh.pPolyScopeSurfaceMesh = pPolyScopeSurfaceMesh;
		}
	}

	void ViewerImpl::registerTetmeshes(const std::vector<TetMeshFEM::SharedPtr>& inTetmeshes)
	{
		tetmeshes = inTetmeshes;
		for (size_t iMesh = 0; iMesh < tetmeshes.size(); iMesh++)
		{
			GAIA::TetMeshFEM::SharedPtr pTetMesh = tetmeshes[iMesh];
			std::string meshName = std::to_string(iMesh);
			MF::STR::padTo(meshName, 4, '0');
			meshName = "Tetmesh_" + meshName;
			polyscope::SurfaceMesh* pPolyScopeSurfaceMesh =
				polyscope::registerSurfaceMesh(meshName, pTetMesh->positions().transpose(), pTetMesh->surfaceFacesTetMeshVIds().transpose());

			polyscopeMeshes.emplace_back();
			PolyScopeSurfaceMesh& polyScopeSurfaceMesh = polyscopeMeshes.back();
			polyScopeSurfaceMesh.meshName = meshName;
			polyScopeSurfaceMesh.meshType = 0;
			polyScopeSurfaceMesh.meshId = iMesh;
			polyScopeSurfaceMesh.pPolyScopeSurfaceMesh = pPolyScopeSurfaceMesh;
		}
	}

}
GAIA::Viewer::Viewer(ViewerParams::SharedPtr pViewerParams_in)
	: pViewerParams(pViewerParams_in)
{
	pImpl = std::make_shared<ViewerImpl>(pViewerParams_in);
}
GAIA::Viewer::~Viewer()
{
}
void GAIA::Viewer::init()
{
	polyscope::options::prepareImGuiFontsCallback = pImpl->setFontcallBack;
	polyscope::options::groundPlaneEnabled = pViewerParams->enableGround;
	polyscope::options::programName = "Gaia Viewer";
	polyscope::options::printPrefix = "[Gaia Viewer]";
	polyscope::init();
	polyscope::state::userCallback = [this]() {
		//std::cout << "Viewer callback invoked." << std::endl;
		this->update(); 
		};
}
void GAIA::Viewer::update()
{
	pImpl->update();
}
void GAIA::Viewer::setAllMeshesToUpdated()
{
	pImpl->setAllMeshesToUpdated();
}
void GAIA::Viewer::show()
{
	std::cout << "Viewer starting." << std::endl;
	polyscope::show();
}
void GAIA::Viewer::registerTrimeshes(const std::vector<TriMeshFEM::SharedPtr>& inTrimeshes)
{
	pImpl->registerTrimeshes(inTrimeshes);
}
void GAIA::Viewer::registerTetmeshes(const std::vector<TetMeshFEM::SharedPtr>& inTetmeshes)
{
	pImpl->registerTetmeshes(inTetmeshes);
}
void GAIA::Viewer::frameTick()
{
	polyscope::frameTick();
}

// !GAIA_NO_GUI
#else
GAIA::Viewer::Viewer(ViewerParams::SharedPtr pViewerParams_in)
{
	GUINoCompliationError();
}
GAIA::Viewer::~Viewer()
{
	GUINoCompliationError();
}
void GAIA::Viewer::frameTick()
{
	GUINoCompliationError();
}
void GAIA::Viewer::update()
{
	GUINoCompliationError();
}
void GAIA::Viewer::registerTrimeshes(const std::vector<TriMeshFEM::SharedPtr>& inTrimeshes)
{
	GUINoCompliationError();
}
void GAIA::Viewer::registerTetmeshes(const std::vector<TetMeshFEM::SharedPtr>& inTetmeshes)
{
	GUINoCompliationError();
}
void GAIA::Viewer::setAllMeshesToUpdated()
{
	GUINoCompliationError();
}
void GAIA::Viewer::show()
{
	GUINoCompliationError();
}
void GAIA::Viewer::init()
{
	GUINoCompliationError();
}
#endif 



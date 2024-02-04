#pragma once
#include "Parser.h"
#include "../VersionTracker/VersionTracker.h"
#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>

template<typename PhysicsFramework>
struct InputHandler {
	template<typename PhysicsCommandParser>
	void handleInput(std::string& inModelInputFile, std::string& inParameterFile, std::string& outFolder,
		PhysicsCommandParser& configs, PhysicsFramework& physics) {

		if (configs.repoRoot != "")
		{
			inModelInputFile = MF::STR::replace(inModelInputFile, "${REPO_ROOT}", configs.repoRoot);
			inParameterFile = MF::STR::replace(inParameterFile, "${REPO_ROOT}", configs.repoRoot);
			outFolder = MF::STR::replace(outFolder, "${REPO_ROOT}", configs.repoRoot);
		}

		if (configs.showGitInfo)
		{
			GitVersionTracker::printGitInfo();
		}

		MF::IO::createFolder(outFolder);

		physics.loadRunningparameters(inModelInputFile, inParameterFile, outFolder);

		if (configs.repoRoot != "")
		{
			for (size_t iObj = 0; iObj < physics.objectParamsList.objectParams.size(); iObj++)
			{
				physics.getObjectParam(iObj).path = MF::STR::replace(physics.getObjectParam(iObj).path, "${REPO_ROOT}",
					configs.repoRoot);
				physics.getObjectParam(iObj).tetsColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).tetsColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);
				physics.getObjectParam(iObj).edgesColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).edgesColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);

			}
			physics.physicsParams().shaderFolderPath = MF::STR::replace(physics.physicsParams().shaderFolderPath, "${REPO_ROOT}",
				configs.repoRoot);

		}
	}
};

template<typename PhysicsFramework>
struct InputHandlerPBDCloth {
	template<typename PhysicsCommandParser>
	void handleInput(std::string& inModelInputFile, std::string& inParameterFile, std::string& outFolder,
		PhysicsCommandParser& configs, PhysicsFramework& physics) {

		if (configs.repoRoot != "")
		{
			inModelInputFile = MF::STR::replace(inModelInputFile, "${REPO_ROOT}", configs.repoRoot);
			inParameterFile = MF::STR::replace(inParameterFile, "${REPO_ROOT}", configs.repoRoot);
			outFolder = MF::STR::replace(outFolder, "${REPO_ROOT}", configs.repoRoot);
		}

		if (configs.showGitInfo)
		{
			GitVersionTracker::printGitInfo();
		}

		MF::IO::createFolder(outFolder);

		physics.loadRunningparameters(inModelInputFile, inParameterFile, outFolder);

		if (configs.repoRoot != "")
		{
			for (size_t iObj = 0; iObj < physics.objectParamsList->size(); iObj++)
			{
				physics.getObjectParam(iObj).path = MF::STR::replace(physics.getObjectParam(iObj).path, "${REPO_ROOT}",
					configs.repoRoot);
				physics.getObjectParam(iObj).triangleColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).tetsColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);
			}
		}
		physics.physicsParams().shaderFolderPath = MF::STR::replace(physics.physicsParams().shaderFolderPath, "${REPO_ROOT}",
			configs.repoRoot);

	}
};

template<typename PhysicsFramework>
struct InputHandlerAPAPCloth {
	template<typename PhysicsCommandParser>
	void handleInput(std::string& inModelInputFile, std::string& inParameterFile, std::string& outFolder,
		PhysicsCommandParser& configs, PhysicsFramework& physics) {

		if (configs.repoRoot != "")
		{
			inModelInputFile = MF::STR::replace(inModelInputFile, "${REPO_ROOT}", configs.repoRoot);
			inParameterFile = MF::STR::replace(inParameterFile, "${REPO_ROOT}", configs.repoRoot);
			outFolder = MF::STR::replace(outFolder, "${REPO_ROOT}", configs.repoRoot);
		}

		if (configs.showGitInfo)
		{
			GitVersionTracker::printGitInfo();
		}

		MF::IO::createFolder(outFolder);

		physics.loadRunningparameters(inModelInputFile, inParameterFile, outFolder);

		if (configs.repoRoot != "")
		{
			for (size_t iObj = 0; iObj < physics.objectParamsList->size(); iObj++)
			{
				physics.getObjectParam(iObj).path = MF::STR::replace(physics.getObjectParam(iObj).path, "${REPO_ROOT}",
					configs.repoRoot);
				physics.getObjectParam(iObj).verticesColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).verticesColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);
			}
		}
		physics.physicsParams().shaderFolderPath = MF::STR::replace(physics.physicsParams().shaderFolderPath, "${REPO_ROOT}",
			configs.repoRoot);

	}
};

template<typename PhysicsFramework>
struct InputHandlerEBD {
	template<typename PhysicsCommandParser>
	void handleInput(std::string& inModelInputFile, std::string& inParameterFile, std::string& outFolder,
		PhysicsCommandParser& configs, PhysicsFramework& physics) {

		if (configs.repoRoot != "")
		{
			inModelInputFile = MF::STR::replace(inModelInputFile, "${REPO_ROOT}", configs.repoRoot);
			inParameterFile = MF::STR::replace(inParameterFile, "${REPO_ROOT}", configs.repoRoot);
			outFolder = MF::STR::replace(outFolder, "${REPO_ROOT}", configs.repoRoot);
		}

		if (configs.showGitInfo)
		{
			GitVersionTracker::printGitInfo();
		}

		MF::IO::createFolder(outFolder);

		physics.loadRunningparameters(inModelInputFile, inParameterFile, outFolder);

		if (configs.repoRoot != "")
		{
			for (size_t iObj = 0; iObj < physics.objectParamsList->objectParams.size(); iObj++)
			{
				physics.getObjectParam(iObj).path = MF::STR::replace(physics.getObjectParam(iObj).path, "${REPO_ROOT}",
					configs.repoRoot);

				physics.getObjectParam(iObj).tetsColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).tetsColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);

				physics.getObjectParam(iObj).verticesColoringCategoriesPath = MF::STR::replace(physics.getObjectParam(iObj).verticesColoringCategoriesPath, "${REPO_ROOT}",
					configs.repoRoot);

				//if (physics.getObjectParam(iObj).materialName == "MassSpring")
				//{
				//	GAIA::ObjectParamsPBDMassSpring::SharedPtr pObjectParamsMaterial = std::static_pointer_cast<GAIA::ObjectParamsPBDMassSpring>(physics.objectParamsList.objectParams[iObj]);
				//	pObjectParamsMaterial->edgesColoringCategoriesPath = MF::STR::replace(pObjectParamsMaterial->edgesColoringCategoriesPath, "${REPO_ROOT}",
				//		configs.repoRoot);

				//}
			}
			physics.physicsParams().shaderFolderPath = MF::STR::replace(physics.physicsParams().shaderFolderPath, "${REPO_ROOT}",
				configs.repoRoot);

		}
	}
};
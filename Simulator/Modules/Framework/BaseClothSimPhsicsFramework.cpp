#pragma once

#include "BaseClothSimPhsicsFramework.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"
#include "../IO/FileIO.h"

using namespace GAIA;

TetMeshFEM::SharedPtr BaseClothPhsicsFramework::initializeMaterial(ObjectParams::SharedPtr objParam, TetMeshMF::SharedPtr pTMeshMF, BasePhysicsParams::SharedPtr physicsParaemters)
{
	std::cout << "Error! initializeMaterial for TetMeshFEM shouldn't be called for cloth simulation!\n";
	std::exit(-1);
	return TetMeshFEM::SharedPtr();
}

void BaseClothPhsicsFramework::initialize()
{
	timeStepEqualize();

	std::cout << "----------------------------------------------------\n"
		<< "Load input mesh files and precomputing topological information.\n"
		<< "----------------------------------------------------\n";

	baseTriMeshes.resize(objectParamsList->objectParams.size(), nullptr);
	//GPUTMeshes.resize(objectParamsList.objectParams.size(), nullptr);
	//#ifdef KEEP_MESHFRAME_MESHES
	//	tMeshesMF.resize(objectParamsList.objectParams.size(), nullptr);
	//#endif // KEEP_MESHFRAME_MESHES

	//cpu_parallel_for(0, objectParamsList.objectParams.size(), [&](int iMesh) {
	for (int iMesh = 0; iMesh < objectParamsList->objectParams.size(); ++iMesh) {


		TriMeshFEM::SharedPtr pTriMesh = initializeMaterial(objectParamsList->objectParams[iMesh], basePhysicsParams, this);
		baseTriMeshes[iMesh] = pTriMesh;

		std::cout << "Added " << objectParamsList->objectParams[iMesh]->materialName
			<< " trimesh: " << objectParamsList->objectParams[iMesh]->path << "\n"
			<< "with " << pTriMesh->numVertices() << " positions and " << pTriMesh->numFaces() << " faces.\n";

	}

	std::cout
		<< "----------------------------------------------------\n"
		<< "Initializing collision detectors. "
		<< "\n----------------------------------------------------" << std::endl;

	// initializeCollisionDetector();
}

void BaseClothPhsicsFramework::writeOutputs(std::string outFolder, int frameId)
{
	std::ostringstream aSs;
	aSs << std::setfill('0') << std::setw(8) << frameId;
	std::string outNumber = aSs.str();
	std::string outFile = outFolder + "/A" + outNumber + "." + basePhysicsParams->outputExt;

	std::string tetMeshOutStatistics = outFolder + "/Statistics";

	if (basePhysicsParams->outputExt == "ply")
	{
		writeAllToPLY(outFile.c_str(), baseTriMeshes, basePhysicsParams->saveAllModelsTogether);
	}
	else if (basePhysicsParams->outputExt == "obj") {
		// writeAllToObj(outFile.c_str(), getSoftBodies(), physicsAllParams.pPhysicsParams.saveAllModelsTogether);
	}
	else
	{
		std::cout << "[Error] Unsupported output EXT name: " << basePhysicsParams->outputExt << std::endl;
		return;
	}

	if (basePhysicsParams->outputRecoveryState && !(frameId % basePhysicsParams->outputRecoveryStateStep))
	{
		std::string stateOutOutPath = outFolder + "/RecoveryStates";
		ClothPhysicsState state;
		std::string outFileState = stateOutOutPath + "/A" + outNumber + ".json";
		state.binary = basePhysicsParams->outputRecoveryStateBinary;
		state.fromPhysics(*this);
		state.writeToJsonFile(outFileState, 2, &outFileState);
	}

}

void GAIA::BaseClothPhsicsFramework::recoverFromState(std::string& stateFile)
{
	ClothPhysicsState state;
	std::cout << "----------------------------------------------------\n"
		<< "Recovering state from:" << stateFile << "\n"
		<< "----------------------------------------------------\n";
	state.loadFromJsonFile(stateFile, &stateFile);
	state.initializePhysics(*this);
	//if (baseCollisionParams->allowDCD && pDCD)
	//{
	//	pDCD->updateBVH(RTC_BUILD_QUALITY_LOW, RTC_BUILD_QUALITY_LOW, true);

	//}
	//if (baseCollisionParams->allowCCD && pCCD)
	//{
	//	pCCD->updateBVH(RTC_BUILD_QUALITY_LOW);
	//}
}

#pragma once

#include "BaseClothSimPhsicsFramework.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"
#include "../IO/FileIO.h"
#include "../SpatialQuery/ColiiderTriMeshSequence.h"

using namespace GAIA;

TetMeshFEM::SharedPtr BaseClothPhsicsFramework::initializeMaterial(ObjectParams::SharedPtr objParam, TetMeshMF::SharedPtr pTMeshMF, BasePhysicsParams::SharedPtr physicsParaemters)
{
	std::cout << "Error! initializeMaterial for TetMeshFEM shouldn't be called for cloth simulation!\n";
	std::exit(-1);
	return TetMeshFEM::SharedPtr();
}

void GAIA::BaseClothPhsicsFramework::parseRunningParameters()
{
	BasePhysicFramework::parseRunningParameters();

	contactDetectorParams = physicsJsonParams["ContactDetectorParams"];
	pClothContactDetectorParameters = std::make_shared<ClothContactDetectorParameters>();
	pClothContactDetectorParameters->fromJson(contactDetectorParams);

	colliderJsonParams = physicsJsonParams["ColliderParams"];
	pDynamicColliderParameter = std::make_shared<DynamicColliderParameters>();
	pDynamicColliderParameter->fromJson(colliderJsonParams);
}

void BaseClothPhsicsFramework::initialize()
{
	timeStepEqualize();

	std::cout << "----------------------------------------------------\n"
		<< "Load input mesh files and precomputing topological information.\n"
		<< "----------------------------------------------------\n";

	baseTriMeshesForSimulation.resize(objectParamsList->objectParams.size(), nullptr);
	//GPUTMeshes.resize(objectParamsList.objectParams.size(), nullptr);
	//#ifdef KEEP_MESHFRAME_MESHES
	//	tMeshesMF.resize(objectParamsList.objectParams.size(), nullptr);
	//#endif // KEEP_MESHFRAME_MESHES

	//cpu_parallel_for(0, objectParamsList.objectParams.size(), [&](int iMesh) {
	for (int iMesh = 0; iMesh < objectParamsList->objectParams.size(); ++iMesh) {


		TriMeshFEM::SharedPtr pTriMesh = initializeMaterial(objectParamsList->objectParams[iMesh], basePhysicsParams, this);
		baseTriMeshesForSimulation[iMesh] = pTriMesh;

		std::cout << "Added " << objectParamsList->objectParams[iMesh]->materialName
			<< " trimesh: " << objectParamsList->objectParams[iMesh]->path << "\n"
			<< "with " << pTriMesh->numVertices() << " positions and " << pTriMesh->numFaces() << " faces.\n";

	}

	std::cout
		<< "----------------------------------------------------\n"
		<< "Initializing collision detectors. "
		<< "\n----------------------------------------------------" << std::endl;

	// initializeCollisionDetector();
	initializeCollider();

	pClothContactDetector = std::make_shared<ClothContactDetector>(pClothContactDetectorParameters);
	triMeshesAll = baseTriMeshesForSimulation;

	for (int i = 0; i < colliderMeshes.size(); i++)
	{
		triMeshesAll.push_back(colliderMeshes[i]);
	}

	pClothContactDetector->initialize(triMeshesAll);
}

void GAIA::BaseClothPhsicsFramework::initializeCollider()
{
	pDynamicCollider = std::make_shared<DynamicCollider>(pDynamicColliderParameter);
	
	nlohmann::json colliderMeshsJson = colliderJsonParams["ColliderMeshes"];
	for (size_t i = 0; i < colliderMeshsJson.size(); i++)
	{
		colliderMeshes.push_back(createColliderMesh(colliderMeshsJson[i]));
	}

	pDynamicCollider->initialize(colliderMeshes);
}



ColliderTrimeshBase::SharedPtr GAIA::BaseClothPhsicsFramework::createColliderMesh(nlohmann::json& colliderMeshJsonParams)
{
	ColliderTrimeshBase::SharedPtr pColliderMesh = nullptr;
	ColliderTrimeshBaseParams::SharedPtr pColliderMeshParams = nullptr;

	if (colliderMeshJsonParams["colliderType"] == "TriMeshSequence")
	{
		pColliderMeshParams = std::make_shared<ColliderTrimeshSequenceParams>();
		pColliderMeshParams->fromJson(colliderMeshJsonParams);
		pColliderMesh = std::make_shared<ColliderTrimeshSequence>();
	}
	else 
	{
		std::cerr << "Error! Unrecognized collider type: " << pColliderMeshParams->colliderType << std::endl;
		exit(-1);
	}

	pColliderMesh->initialize(pColliderMeshParams);
	return pColliderMesh;
}

void GAIA::BaseClothPhsicsFramework::initializeViewer()
{
	BasePhysicFramework::initializeViewer();
	pViewer->registerTrimeshes(triMeshesAll);

}

void GAIA::BaseClothPhsicsFramework::updateCollider()
{
	pDynamicCollider->updateColliderMeshes(frameId, substep, iIter, basePhysicsParams->numSubsteps, basePhysicsParams->iterations);
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
		writeAllToPLY(outFile.c_str(), baseTriMeshesForSimulation, basePhysicsParams->saveAllModelsTogether);
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

void GAIA::BaseClothPhsicsFramework::setSimulatedMeshToUpToDateStatus(bool updated)
{
	for (IdType iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		if (baseTriMeshesForSimulation[iMesh]->activeForSim) {
			baseTriMeshesForSimulation[iMesh]->updated = updated;
		}
	}
}




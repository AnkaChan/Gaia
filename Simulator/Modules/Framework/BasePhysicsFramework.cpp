#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>
#include <MeshFrame/Utility/Time.h>
#include "BasePhysicsFramework.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

#include "../CollisionDetector/DiscreteCollisionDetector.h"
#include "../CollisionDetector/ContinuousCollisionDetector.h"


#include "../IO/FileIO.h"
#include "../Timer/Timer.h"

#include "../Parallelization/CPUParallelization.h"

#include "../Viewer/Viewer.h"

using namespace GAIA;

size_t GAIA::ObjectParamsList::size()
{
	return objectParams.size();
}

ObjectParams& GAIA::ObjectParamsList::getObjectParam(int iObj)
{
	return *objectParams[iObj];
}

bool GAIA::ObjectParamsList::fromJson(nlohmann::json& objectParam)
{
	nlohmann::json modelsInfo = objectParam["Models"];

	for (nlohmann::json& model : modelsInfo) {
		std::string materialName;
		parseJsonParameters(model, "materialName", materialName);

		ObjectParams::SharedPtr pObjParams = createObjectParam(materialName);

		if (materialName == "NeoHookean")
		{
			pObjParams->materialType = NeoHookean;

		}
		else if (materialName == "MassSpring") {
			pObjParams->materialType = MassSpring;

		}
		else if (materialName == "StVK_triMesh")
		{
			pObjParams->materialType = StVK_triMesh;
		}


		pObjParams->fromJson(model);
		pObjParams->materialName = materialName;
		objectParams.push_back(pObjParams);
	}
	return true;
}

bool GAIA::ObjectParamsList::toJson(nlohmann::json& objectParam)
{
	for (size_t iObj = 0; iObj < objectParams.size(); iObj++)
	{
		nlohmann::json obj;
		objectParams[iObj]->toJson(obj);
		objectParam["Models"].push_back(obj);
	}
	return true;
}

void GAIA::BasePhysicFramework::updateWorldBox()
{
	for (int iMesh = 0; iMesh < basetetMeshes.size(); iMesh++)
	{
		bool changed = false;
		TetMeshFEM* pTM = basetetMeshes[iMesh].get();
		for (int iDim = 0; iDim < 3; iDim++)
		{
			FloatingType dimMin = pTM->positions().row(iDim).minCoeff();
			FloatingType dimMax = pTM->positions().row(iDim).maxCoeff();
			if (basePhysicsParams->worldBounds(iDim, 0) > dimMin)
			{
				basePhysicsParams->worldBounds(iDim, 0) = dimMin;
				changed = true;
			}

			if (basePhysicsParams->worldBounds(iDim, 1) < dimMax)
			{
				basePhysicsParams->worldBounds(iDim, 1) = dimMax;
				changed = true;
			}
		}

		std::cout << "World box updates to: \n" << basePhysicsParams->worldBounds.transpose() << std::endl;
	}
}


void GAIA::BasePhysicFramework::initialize()
{
	timeStepEqualize();

	std::cout << "----------------------------------------------------\n"
		<< "Load input mesh files and precomputing topological information.\n"
		<< "----------------------------------------------------\n";
	basetetMeshes.resize(objectParamsList->objectParams.size(), nullptr);
	//GPUTMeshes.resize(objectParamsList.objectParams.size(), nullptr);
#ifdef KEEP_MESHFRAME_MESHES
	tMeshesMF.resize(objectParamsList.objectParams.size(), nullptr);
#endif // KEEP_MESHFRAME_MESHES

	std::map<std::string, TetMeshMF::SharedPtr> tmeshMFPtrs;
	std::mutex tmeshMFPtrs_lock;

	numAllVertices = 0;
	numAllTets = 0;
	numAllEdges = 0;

	cpu_parallel_for(0, objectParamsList->objectParams.size(), [&](int iMesh) {
		//for (int iMesh = 0; iMesh < objectParamsList->objectParams.size(); ++iMesh) {
		TetMeshMF::SharedPtr pTM_MF;
		std::string modelPath = objectParamsList->objectParams[iMesh]->path;
		MF::IO::FileParts fp = MF::IO::fileparts(modelPath);

		// scope of meshMFPtrsLockGuard

		bool loadSucceed = false;
		{
			std::lock_guard<std::mutex> meshMFPtrsLockGuard(tmeshMFPtrs_lock);

			auto pTMeshMFItem = tmeshMFPtrs.find(modelPath);
			if (pTMeshMFItem == tmeshMFPtrs.end())
			{
				pTM_MF = std::make_shared<TetMeshMF>();
				if (fp.ext == ".t")
				{
					//pTM->_load_t(param.inputModelPath.c_str(), true);
					pTM_MF->load_t(objectParamsList->objectParams[iMesh]->path.c_str());
					loadSucceed = true;
					tmeshMFPtrs.insert({ modelPath, pTM_MF });
				}
				else if (fp.ext == ".vtk") {
					std::cout << "Currently vtk file is not supported! " << std::endl;
					loadSucceed = false;
					//continue;
				}
				else
				{
					std::cout << "Unsupported file format: " << fp.ext << std::endl;
					loadSucceed = false;
					//continue;
				}
			}
			else
			{
				loadSucceed = true;
				pTM_MF = pTMeshMFItem->second;
			}

			std::cout << "Adding " << objectParamsList->objectParams[iMesh]->materialName
				<< " tetmesh: " << modelPath << "\n"
				<< "with " << pTM_MF->numVertices() << " vertices and " << pTM_MF->numTets() << " tets.\n";
		}

		if (loadSucceed) {
			TetMeshFEM::SharedPtr pTMesh = initializeMaterial(objectParamsList->objectParams[iMesh], pTM_MF, basePhysicsParams);
			pTMesh->meshId = iMesh;
			basetetMeshes[iMesh] = pTMesh;

			numAllVertices += pTMesh->numVertices();
			numAllTets += pTMesh->numTets();
			numAllEdges += pTMesh->numEdges();

		}
		else
		{
			std::cout << "[Error] Fail to load: " << objectParamsList->objectParams[iMesh]->path << std::endl;
			exit(-1);
		}
		});

	std::cout
		<< "----------------------------------------------------\n"
		<< "Initializing collision detectors. "
		<< "\n----------------------------------------------------" << std::endl;

	initializeCollisionDetector();

	if (basePhysicsParams->checkAndUpdateWorldBounds)
	{
		updateWorldBox();
	}

	disableModelsLatterToAppear();

	if (pViewerParams->enableViewer)
	{
		initializeViewer();
	}
}

void GAIA::BasePhysicFramework::initializeViewer()
{
	pViewer = std::make_shared<Viewer>(pViewerParams);
	pViewer->init();
	pViewer->registerTetmeshes(basetetMeshes);
}

void GAIA::BasePhysicFramework::disableModelsLatterToAppear()
{
	for (int iMesh = 0; iMesh < numMeshes(); iMesh++)
	{
		TetMeshFEM* pTetMesh = basetetMeshes[iMesh].get();

		if (pTetMesh->pObjectParams->frameToAppear > 0)
		{
			pTetMesh->activeForCollision = false;
			pTetMesh->activeForMaterialSolve = false;
		}
	}
}

void GAIA::BasePhysicFramework::enableModels()
{
	for (int iMesh = 0; iMesh < numMeshes(); iMesh++)
	{
		TetMeshFEM* pTetMesh = basetetMeshes[iMesh].get();

		if (pTetMesh->pObjectParams->frameToAppear >= frameId) {

			if (false == pTetMesh->activeForCollision)
			{
				pTetMesh->activeForCollision = true;
			}
			if (false == pTetMesh->activeForMaterialSolve)
			{
				pTetMesh->activeForMaterialSolve = true;
			}
		}
	}
}

void GAIA::BasePhysicFramework::timeStepEqualize()
{
	basePhysicsParams->dt = basePhysicsParams->timeStep / basePhysicsParams->numSubsteps;
	basePhysicsParams->dtSqrReciprocal = 1.f / (basePhysicsParams->dt * basePhysicsParams->dt);
}

void GAIA::BasePhysicFramework::initializeCollisionDetector()
{
	if (baseCollisionParams->allowCCD)
	{
		pCCD = std::make_shared<ContinuousCollisionDetector>(*baseCollisionParams);
		pCCD->initialize(basetetMeshes);
	}
	if (baseCollisionParams->allowDCD) {
		pDCD = std::make_shared<DiscreteCollisionDetector>(*baseCollisionParams);
		pDCD->initialize(basetetMeshes);
	}
}

void GAIA::BasePhysicFramework::recoverFromState(std::string& stateFile)
{
	PhysicsState state;
	std::cout << "----------------------------------------------------\n"
		<< "Recovering state from:" << stateFile << "\n"
		<< "----------------------------------------------------\n";
	state.loadFromJsonFile(stateFile, &stateFile);
	state.initializePhysics(*this);

	if (baseCollisionParams->allowDCD)
	{
		pDCD->updateBVH(RTC_BUILD_QUALITY_LOW, RTC_BUILD_QUALITY_LOW, true);

	}
	if (baseCollisionParams->allowCCD)
	{
		pCCD->updateBVH(RTC_BUILD_QUALITY_LOW);
	}
}

void GAIA::BasePhysicFramework::simulate()
{

	setUpOutputFolders(outputFolder);
	std::cout
		<< "----------------------------------------------------\n"
		<< "Output folder is: " << outputFolder << std::endl;

	writeOutputs(outputFolder, frameId);
	std::cout
		<< "----------------------------------------------------\n"
		<< "Starting Sims\n"
		<< "----------------------------------------------------"
		<< std::endl;

	baseTimeStatistics->setToZero();

	while (frameId < basePhysicsParams->numFrames) {
		TICK(timeCsmpFrame);
		debugOperation(DEBUG_LVL_INFO, [&]() {
			std::cout
				<< "----------------------------------------------------\n"
				<< "Frame " << frameId + 1 << " begin.\n"
				;
			});

		enableModels();
		runStep();

		++frameId;

		if (basePhysicsParams->saveOutputs)
		{
			TICK(timeCsmpSaveOutputs);
			writeOutputs(outputFolder, frameId + 1);
			TOCK_STRUCT((*baseTimeStatistics), timeCsmpSaveOutputs);
		}

		if (pViewerParams->enableViewer)
		{
			pViewer->update();
			pViewer->frameTick();
		}

		TOCK_STRUCT((*baseTimeStatistics), timeCsmpFrame);

		debugPrint(DEBUG_LVL_INFO, baseTimeStatistics->getString());
		debugOperation(DEBUG_LVL_INFO, [&]() {
			std::cout
				<< "Frame " << frameId + 1 << " completed, Time consumption: " << baseTimeStatistics->timeCsmpFrame << "\n"
				<< "----------------------------------------------------\n";
			});

		baseTimeStatistics->setToZero();
	}
}

std::string GAIA::BasePhysicFramework::getDebugFolder()
{
	if (basePhysicsParams->debugOutFolder != "")
	{
		return basePhysicsParams->debugOutFolder;
	}
	else
	{
		std::string outputDbgFolder = outputFolder + "/Debug";
		MF::IO::createFolder(outputDbgFolder);
		return outputDbgFolder;
	}
}

std::shared_ptr<CollisionDetectionParamters> GAIA::BasePhysicFramework::createCollisionParams()
{
	return std::make_shared<CollisionDetectionParamters>();
}

void GAIA::BasePhysicFramework::setUpOutputFolders(std::string outFolder)
{
	MF::IO::createFolder(outFolder);

	std::string tetMeshOutOutPath = outFolder + "/TetMesh";
	if (basePhysicsParams->outputVTK || basePhysicsParams->outputT)
	{
		MF::IO::createFolder(tetMeshOutOutPath);
	}

	if (basePhysicsParams->saveSimulationParameters)
	{
		std::string paramsOutOutPath = outFolder + "/SimulationParameters_" + MF::Time::getCurTimeString();
		MF::IO::createFolder(paramsOutOutPath);
		saveExperimentParameters(paramsOutOutPath);
	}

	if (basePhysicsParams->outputStatistics) {
		std::string statisticsOutOutPath = outFolder + "/Statistics";
		MF::IO::createFolder(statisticsOutOutPath);
	}

	if (basePhysicsParams->outputRecoveryState) {
		std::string stateOutOutPath = outFolder + "/RecoveryStates";
		MF::IO::createFolder(stateOutOutPath);
	}
}

void GAIA::BasePhysicFramework::writeOutputs(std::string outFolder, int frameId)
{
	std::ostringstream aSs;
	aSs << std::setfill('0') << std::setw(8) << frameId;
	std::string outNumber = aSs.str();
	std::string outFile = outFolder + "/A" + outNumber + "." + basePhysicsParams->outputExt;

	std::string tetMeshOutOutPath = outFolder + "/TetMesh";
	std::string tetMeshOutStatistics = outFolder + "/Statistics";

	if (basePhysicsParams->outputExt == "bin")
	{
		if (frameId == 0)
			// 
			//  template
		{
			std::string templateMeshOutOutPath = outFolder + "/TemplateMesh";
			MF::IO::createFolder(templateMeshOutOutPath);

			std::string templateMeshesName = outFolder + "/TemplateMesh/TemplateMesh.ply";
			writeAllToPLY(templateMeshesName.c_str(), basetetMeshes, false, true);
		}
		writeAllToBinary(outFile.c_str(), basetetMeshes);

		if (!(frameId % basePhysicsParams->binaryModeVisualizationSteps))
		{
			std::string outFileVis = outFolder + "/A" + outNumber + ".ply";

			writeAllToPLY(outFileVis.c_str(), basetetMeshes, basePhysicsParams->saveAllModelsTogether, basePhysicsParams->saveAllModelsTogether);
		}
	}
	else if (basePhysicsParams->outputExt == "ply")
	{
		writeAllToPLY(outFile.c_str(), basetetMeshes, basePhysicsParams->saveAllModelsTogether, basePhysicsParams->saveAllModelsTogether);
	}
	else if (basePhysicsParams->outputExt == "obj") {
		// writeAllToObj(outFile.c_str(), getSoftBodies(), physicsAllParams.pPhysicsParams.saveAllModelsTogether);
	}
	else
	{
		std::cout << "[Error] Unsupported output EXT name: " << basePhysicsParams->outputExt << std::endl;
		return;
	}

	if (basePhysicsParams->outputT) {
		std::string outFileT = tetMeshOutOutPath + "/A" + outNumber + ".t";
		//writeAllTetMeshToT(outFileT.c_str(), getSoftBodies());
	}

	if (basePhysicsParams->outputVTK) {
		std::string outFileVtk = tetMeshOutOutPath + "/A" + outNumber + ".vtk";
		// writeAllTetMeshToVtk(outFileVtk.c_str(), getSoftBodies());
	}

	if (basePhysicsParams->outputStatistics)
	{
		//pSoftbodyManager->statistics.writeToJsonFile(outFolder + "/Statistics.json");
		baseTimeStatistics->writeToJsonFile(tetMeshOutStatistics + "/A" + outNumber + ".json");
		if (frameId % basePhysicsParams->outputRecoveryStateStep == 0) {
#ifdef DO_COLLISION_STATISTICS
			collisionStatistics.writeToJsonFile(tetMeshOutStatistics + "/CollisionStatisticsAll.json");
#endif // DO_COLLISION_STATISTICS

		}
	}

	if (basePhysicsParams->outputRecoveryState && !(frameId % basePhysicsParams->outputRecoveryStateStep))
	{
		std::string stateOutOutPath = outFolder + "/RecoveryStates";
		PhysicsState state;

		std::string outFileState = stateOutOutPath + "/A" + outNumber + ".json";
		state.binary = basePhysicsParams->outputRecoveryStateBinary;
		state.fromPhysics(*this);
		state.writeToJsonFile(outFileState, 2, &outFileState);

	}
}

bool GAIA::BasePhysicFramework::writeSimulationParameters(nlohmann::json& outPhysicsParams)
{
	basePhysicsParams->toJson(outPhysicsParams["PhysicsParams"]);
	baseCollisionParams->toJson(outPhysicsParams["CollisionParams"]);
	pViewerParams->toJson(physicsJsonParams["ViewerParams"]);
	return true;

}

void GAIA::BasePhysicFramework::debugOperation(int debugLvl, std::function<void()> ops)
{
	debugInfoGen(basePhysicsParams->debugVerboseLvl, debugLvl, ops);
}

void GAIA::BasePhysicFramework::debugPrint(int debugLvl, std::string info)
{
	debugOperation(debugLvl, [&]() {
		std::cout << info;
		});
}

void GAIA::BasePhysicFramework::saveDebugState(const std::string customName, bool saveMesh, const std::string outfolder)
{
	std::string outName = outfolder;
	if (outName == "") 
	{
		outName = getDebugFolder();
	}
	std::ostringstream aSs;
	aSs << customName << "Frame" << std::setfill('0') << std::setw(8) << frameId << "_iter" << std::setw(4) << iIter;
	std::string outNumber = aSs.str();
	std::string outNameJson = outName + "/" + outNumber + ".json";


	PhysicsState state;

	state.binary = basePhysicsParams->outputRecoveryStateBinary;
	state.fromPhysics(*this);
	state.writeToJsonFile(outNameJson, 2, &outNameJson);

	if (saveMesh)
	{
		std::string outNameMesh = outName + "/" + outNumber + ".ply";
		writeAllToPLY(outNameMesh.c_str(), basetetMeshes, basePhysicsParams->saveAllModelsTogether, false);
	}
}


void BasePhysicFramework::saveExperimentParameters(const std::string& paramsOutOutPath, int indent)
{
	objectParamsList->writeToJsonFile(paramsOutOutPath + "/Models.json", 2);

	//physicsAllParams.writeToJsonFile();
	nlohmann::json outPhysicsParams;
	std::ofstream ofs(paramsOutOutPath + "/Parameters.json");
	if (ofs.is_open())
	{
		writeSimulationParameters(outPhysicsParams);
		ofs << outPhysicsParams.dump(2);
	}
	else
	{
		std::cout << "[Error] Fail to open file: " << paramsOutOutPath << std::endl;
	}

}


void BasePhysicFramework::loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder)
{
	std::cout << "----------------------------------------------------\n"
		<< "Loading input parameter files.\n"
		<< "----------------------------------------------------\n";

	outputFolder = outFolder;
	inputModelListFile = inModelInputFile;
	inputParamFile = inParameterFile;

	MF::loadJson(inModelInputFile, modelJsonParams);
	MF::loadJson(inParameterFile, physicsJsonParams);

	parseRunningParameters(modelJsonParams, physicsJsonParams);

}

void GAIA::BasePhysicFramework::parseRunningParameters(nlohmann::json& inModelParams, nlohmann::json& inPhysicsParams)
{
	basePhysicsParams = createPhysicsParams();
	baseCollisionParams = createCollisionParams();

	basePhysicsParams->fromJson(inPhysicsParams["PhysicsParams"]);
	baseCollisionParams->fromJson(inPhysicsParams["CollisionParams"]);

	objectParamsList = createObjectParamsList();
	objectParamsList->fromJson(inModelParams);

	baseTimeStatistics = createRunningTimeStatistics();

	pViewerParams = std::make_shared<ViewerParams>();
	pViewerParams->fromJson(physicsJsonParams["ViewerParams"]);
}
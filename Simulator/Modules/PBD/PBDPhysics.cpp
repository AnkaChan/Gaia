#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>
#include <MeshFrame/Utility/Time.h>
#include "../CollisionDetector/DiscreteCollisionDetector.h"
#include "../CollisionDetector/ContinuousCollisionDetector.h"
#include "../CollisionDetector/VolumetricCollisionDetector.h"

#include <memory>
#include "PBDPhysics.h"

#include "PBDTetMeshGeneralCompute.h"

#include "../Parallelization/CPUParallelization.h"

#include "../Timer/Timer.h"
#include "../Timer/RunningTimeStatistics.h"

#include "../IO/FileIO.h"

#include <unordered_set>


// #include "../CollisionDetector/DiscreteCollisionDetector.h"

void GAIA::PBDPhysics::loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder)
{

	std::cout << "----------------------------------------------------\n"
		<< "Loading input parameter files.\n"
		<< "----------------------------------------------------\n";

	physicsAllParams.loadFromJsonFile(inParameterFile);
	objectParamsList.loadFromJsonFile(inModelInputFile);


	inputModelListFile = inModelInputFile;
	inputParamFile = inParameterFile;
	outputFolder = outFolder;

	dt = physicsAllParams.physicsParams.timeStep / physicsAllParams.physicsParams.numSubsteps;
	
	// make parameters step invariant
	for (size_t iMesh = 0; iMesh < objectParamsList.objectParams.size(); iMesh++)
	{
		objectParamsList.getObjectParam(iMesh).exponentialVelDamping = std::pow(objectParamsList.getObjectParam(iMesh).exponentialVelDamping, dt);
		objectParamsList.getObjectParam(iMesh).constantVelDamping = objectParamsList.getObjectParam(iMesh).constantVelDamping * dt;

	}

}

void GAIA::PBDPhysics::saveExperimentParameters(const std::string& paramsOutOutPath)
{
	objectParamsList.writeToJsonFile(paramsOutOutPath + "/Models.json", 2);
	physicsAllParams.writeToJsonFile(paramsOutOutPath + "/Parameters.json", 2);
}

void GAIA::PBDPhysics::syncAllToGPU(bool sync)
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->syncToGPU(false, cudaStreams[iMesh]);
	}
	if (sync)
	{
		cudaDeviceSynchronize();
	}
}

void GAIA::PBDPhysics::syncAllToCPU(bool sync)
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->syncToCPU(false, cudaStreams[iMesh]);
	}
	if (sync)
	{
		cudaDeviceSynchronize();
	}
}

void GAIA::PBDPhysics::syncAllToCPUVertPosOnly(bool sync)
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->syncToCPUVertPosOnly(false, cudaStreams[iMesh]);
	}
	if (sync)
	{
		cudaDeviceSynchronize();
	}
}

void GAIA::PBDPhysics::syncAllToInvertedSignOnly(bool sync)
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->syncToCPUInversionSignOnly(false, cudaStreams[iMesh]);
	}
	if (sync)
	{
		cudaDeviceSynchronize();
	}

	//PBDTetMeshFEM* pTM = tMeshes[0].get();
	//int numInverted = 0;
	//for (int iTet = 0; iTet < pTM->numTets(); iTet++)
	//{
	//	if (pTM->tetsInvertedSign[iTet]) {
	//		++numInverted;
	//	}
	//}
	//std::cout << "num inverted after solving: " << numInverted << std::endl;

	//PBDTetMeshFEM* pTM = tMeshes[0].get();
	//int numInvertedVerts = 0;
	//for (int iV = 0; iV < pTM->numVertices(); iV++)
	//{
	//	if (pTM->verticesInvertedSign[iV]) {
	//		++numInvertedVerts;
	//	}
	//}
	//std::cout << "num inverted verts after solving: " << numInvertedVerts << std::endl;
}


void GAIA::PBDPhysics::syncAllToGPUVertPosOnly(bool sync)
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->syncToGPUVertPosOnly(false, cudaStreams[iMesh]);
	}
	if (sync)
	{
		cudaDeviceSynchronize();
	}
}

void GAIA::PBDPhysics::recoverFromState(std::string& stateFile)
{
	PBDPhysicsState state;
	std::cout << "----------------------------------------------------\n"
		<< "Recovering state from:" << stateFile << "\n"
		<< "----------------------------------------------------\n";
	state.loadFromJsonFile(stateFile);
	state.initializeSoftBodyManager(*this, frameId);

	pDCD->updateBVH(RTC_BUILD_QUALITY_LOW, RTC_BUILD_QUALITY_LOW, true);
	pCCD->updateBVH(RTC_BUILD_QUALITY_LOW);
}

bool GAIA::PBDPhysics::initializeGPU()
{
	std::cout << "----------------------------------------------------\n"
		<< "Load input mesh files and precomputing topological information.\n"
		<< "----------------------------------------------------\n";
	tMeshes.resize(objectParamsList.objectParams.size(), nullptr);
	GPUTMeshes.resize(objectParamsList.objectParams.size(), nullptr);
#ifdef KEEP_MESHFRAME_MESHES
	tMeshesMF.resize(objectParamsList.objectParams.size(), nullptr);
#endif // KEEP_MESHFRAME_MESHES

	//cpu_parallel_for(0, objectParamsList.objectParams.size(), [&](int iMesh) {
	for (int iMesh = 0; iMesh < objectParamsList.objectParams.size(); ++iMesh) {
		TetMeshMF::SharedPtr pTM_MF = std::make_shared<TetMeshMF>();
		MF::IO::FileParts fp = MF::IO::fileparts(objectParamsList.objectParams[iMesh]->path);

		bool loadSucceed = false;
		if (fp.ext == ".t")
		{
			//pTM->_load_t(param.inputModelPath.c_str(), true);
			pTM_MF->load_t(objectParamsList.objectParams[iMesh]->path.c_str());
			loadSucceed = true;
		}
		else if (fp.ext == ".vtk") {
			std::cout << "Currently vtk file is not supported! " << std::endl;
			//continue;
		}
		else
		{
			std::cout << "Unsupported file format: " << fp.ext << std::endl;
			//continue;
		}

		if (loadSucceed) {
#ifdef KEEP_MESHFRAME_MESHES
				tMeshesMF[iMesh] = pTM_MF;
#endif // KEEP_MESHFRAME_MESHES

			switch (objectParamsList.objectParams[iMesh]->materialType)
			{
			case NeoHookean:
			{
				PBDTetMeshNeoHookean::SharedPtr pTetMeshNeoHookean = std::make_shared<PBDTetMeshNeoHookean>();
				tMeshes[iMesh] = pTetMeshNeoHookean;
				pTetMeshNeoHookean->initialize(objectParamsList.objectParams[iMesh], pTM_MF, this);
				GPUTMeshes[iMesh] = pTetMeshNeoHookean->getTetMeshGPU();
				break;
			}
			case MassSpring:
			{
				PBDTetMeshMassSpring::SharedPtr pTetMeshMassSpring = std::make_shared<PBDTetMeshMassSpring>();
				tMeshes[iMesh] = pTetMeshMassSpring;
				pTetMeshMassSpring->initialize(objectParamsList.objectParams[iMesh], pTM_MF, this);
				GPUTMeshes[iMesh] = pTetMeshMassSpring->getTetMeshGPU();
			}
				break;
			default:
				break;
			}
		}
		else
		{
			std::cout << "[Error] Fail to load: " << objectParamsList.objectParams[iMesh]->path << std::endl;
			exit(-1);
		}
	}


	std::cout
		<< "----------------------------------------------------\n"
		<< "Initializing collision detectors. "
		<< "\n----------------------------------------------------" << std::endl;

	pDCD = std::make_shared<DiscreteCollisionDetector>(physicsAllParams.collisionParams);
	pCCD = std::make_shared<ContinuousCollisionDetector>(physicsAllParams.collisionParams);

	ccdResultsAll.resize(tMeshes.size());
	dcdResultsAll.resize(tMeshes.size());

	initialPenetratedVertIds.resize(tMeshes.size());
	initialPenetrationFreeVertIds.resize(tMeshes.size());

	std::vector<std::shared_ptr<TetMeshFEM>> tMeshPtrsBase;
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshPtrsBase.push_back(tMeshes[iMesh]);
		ccdResultsAll[iMesh].resize(tMeshes[iMesh]->numSurfaceVerts());
		dcdResultsAll[iMesh].resize(tMeshes[iMesh]->numSurfaceVerts());
	}
	pDCD->initialize(tMeshPtrsBase);
	pCCD->initialize(tMeshPtrsBase);

	if (physicsAllParams.collisionParams.allowVolumetricCollision)
	{
		pVolCD = std::make_shared<VolumetricCollisionDetector>(physicsAllParams.collisionParams.volCollisionParams);
		pVolCD->initialize(tMeshPtrsBase);
	}

	// initialize GPUS
	cudaStreams.resize(tMeshes.size());
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		cudaStreamCreate(&cudaStreams[iMesh]);

	}
	return true;

}

void GAIA::PBDPhysics::updateGPUMeshes()
{
	for (int iMesh = 0; iMesh < tMeshes.size(); ++iMesh)
	{
		GAIA::PBDTetMeshFEM* pTM = tMeshes[iMesh].get();

		GAIA::PBDTetMeshNeoHookean* pTMNeoHookean;
		GAIA::PBDTetMeshMassSpring* pTMMassSpring;
		switch (pTM->pObjectParams->materialType)
		{
		case GAIA::NeoHookean:
			pTMNeoHookean = (GAIA::PBDTetMeshNeoHookean*)pTM;
			pTMNeoHookean->copyAllToGPU();
			pTMNeoHookean->pTetMeshGPUBuffer->fromCPU(&pTMNeoHookean->tetMeshGPU);

		default:
			pTMMassSpring = (GAIA::PBDTetMeshMassSpring*)pTM;
			pTMMassSpring->copyAllToGPU();
			pTMMassSpring->pTetMeshGPUBuffer->fromCPU(&pTMMassSpring->tetMeshGPU);
			break;
		}
	}
}

void GAIA::PBDPhysics::simulateGPU()
{
	if (physicsParams().checkAndUpdateWorldBounds)
	{
		updateWorldBox();
	}
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

	timeStatistics.setToZero();

	while (frameId < physicsAllParams.physicsParams.numFrames) {
		TICK(timeCsmpFrame);
		runStepGPU();

		TICK(timeCsmpSaveOutputs);
		writeOutputs(outputFolder, frameId + 1);
		TOCK_STRUCT(timeStatistics, timeCsmpSaveOutputs);

		TOCK_STRUCT(timeStatistics, timeCsmpFrame);
		
		std::cout
			<< "----------------------------------------------------\n"
			<< "Frame " << frameId + 1 << " completed, Time consumption: " << timeStatistics.timeCsmpFrame
			<< "\n----------------------------------------------------" << std::endl;
		timeStatistics.print();

		++frameId;
	}


}

void GAIA::PBDPhysics::simulateGPU_debugOnCPU()
{
	if (physicsParams().checkAndUpdateWorldBounds)
	{
		updateWorldBox();
	}

	std::vector<PBDTetMeshFEMGPU*> gpuTMeshesOnCPU;
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		switch (objectParamsList.objectParams[iMesh]->materialType)
		{
		case NeoHookean:
		{
			gpuTMeshesOnCPU.push_back(&(std::static_pointer_cast<PBDTetMeshNeoHookean>(tMeshes[iMesh])->tetMeshGPU_forCPU));
			break;
		}
		case MassSpring:
			gpuTMeshesOnCPU.push_back(&(std::static_pointer_cast<PBDTetMeshMassSpring>(tMeshes[iMesh])->tetMeshGPU_forCPU));
			break;
		default:
			break;
		}
	}

	setUpOutputFolders(outputFolder);
	writeOutputs(outputFolder, frameId);
	std::cout
		<< "----------------------------------------------------\n"
		<< "Starting Sims\n"
		<< "----------------------------------------------------"
		<< std::endl;

	timeStatistics.setToZero();

	while (frameId < physicsAllParams.physicsParams.numFrames) {
		TICK(timeCsmpFrame);
		runStepGPU_debugOnCPU(gpuTMeshesOnCPU);

		TICK(timeCsmpSaveOutputs);
		writeOutputs(outputFolder, frameId + 1);
		TOCK_STRUCT(timeStatistics, timeCsmpSaveOutputs);

		TOCK_STRUCT(timeStatistics, timeCsmpFrame);

		std::cout
			<< "----------------------------------------------------\n"
			<< "Frame " << frameId + 1 << " completed, Time consumption: " << timeStatistics.timeCsmpFrame
			<< "\n----------------------------------------------------" << std::endl;
		timeStatistics.print();

		++frameId;
	}
}

void GAIA::PBDPhysics::runStepGPU()
{	
	timeStatistics.setToZero();

	collisionDetectionCounter = 0;

	for (substep = 0; substep < physicsParams().numSubsteps; substep++) {
		if (physicsParams().showSubstepProgress) {
			std::cout << "---Substep #" << substep << " | time consumption till now: "
				<< timeStatistics.timeCsmpAllSubSteps << std::endl;
		}
		TICK(timeCsmpAllSubSteps);

		TICK(timeCsmpInitialCollision);
		// do the initial detection to categorize the vertices into 2 types: penetration free ones and penetrated ones

		if (!(substep % physicsParams().collisionDetectionSubSteps))
		{
			// no need to build BVH here, it will be rebuilt in the detection part
			initialCollisionsDetection();
		}

		TOCK_STRUCT(timeStatistics, timeCsmpInitialCollision);
		// run on CPU

		// std::cout << tMeshes[0]->vertPos().transpose() << std::endl;

		// run on GPU
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			applyDeformers();

			TICK(timeCsmpMaterialSolve);
			if (iIter != 0) {
				syncAllToGPUVertPosOnly(false);
			}
			else
			{
				applyInitialGuesses();
				evaluateConvergence(substep, 0);
				syncAllToGPU(false);
			}

			materialBoundarySolveGPU();

			TICK(timeCsmpInversionSolve);
			InversionSolveGPU();
			TOCK_STRUCT(timeStatistics, timeCsmpInversionSolve);

			syncAllToCPUVertPosOnly(false);
			syncAllToInvertedSignOnly(true);

			evaluateConvergence(substep, iIter+1);

			//cudaDeviceSynchronize();

			TOCK_STRUCT(timeStatistics, timeCsmpMaterialSolve);

			collisionSolve();
			
			applyConstraints();

			applyPostSolvingDeformers();

			//cudaDeviceSynchronize();
		}
		//evaluateActiveness();

		// collision solve
		// run on CPU
		//syncAllToCPU();

		// std::cout << tMeshes[0]->vertPos().transpose() << std::endl;

		updateVelocities();
		++collisionDetectionCounter;
		curTime += dt;
		TOCK_STRUCT(timeStatistics, timeCsmpAllSubSteps);
	}
}

void GAIA::PBDPhysics::runStepGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU)
{
	timeStatistics.setToZero();
	TICK(timeCsmpAllSubSteps);
	int collisionDetectionCounter = 0;

	for (substep = 0; substep < physicsParams().numSubsteps; substep++) {
		if (physicsParams().showSubstepProgress) {
			std::cout << "---Substep #" << substep << std::endl;
		}

		TICK(timeCsmpInitialCollision);
		// do the initial detection to categorize the vertices into 2 types: penetration free ones and penetrated ones

		if (!(substep % physicsParams().collisionDetectionSubSteps))
		{
			// no need to build BVH here, it will be rebuilt in the detection part
			initialCollisionsDetection();
		}
		applyInitialGuesses();
		TOCK_STRUCT(timeStatistics, timeCsmpInitialCollision);
		// run on CPU

		// syncAllToGPU();
		// std::cout << tMeshes[0]->vertPos().transpose() << std::endl;

		// run GPU code on CPU
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			TICK(timeCsmpMaterialSolve);
			materialBoundarySolveGPU_debugOnCPU(gpuTMeshesOnCPU);
			TOCK_STRUCT(timeStatistics, timeCsmpMaterialSolve);

			TICK(timeCsmpInversionSolve);
			solveInvertedTets();
			TOCK_STRUCT(timeStatistics, timeCsmpInversionSolve);

			TICK(timeCsmpUpdatingCollisionInfoDCD);
			if((collisionDetectionCounter >= physicsParams().collisionDetectionSubSteps || collisionDetectionCounter == 0)
				&& ((!physicsParams().doCollDetectionOnlyForFirstIteration) || (iIter == 0)))
			{
				if (collisionDetectionCounter != 0)
				{
					collisionDetectionCounter = 0;
				}

				positiveCollisionDetectionResults.clear();

				bool rebuildCCDBVH = (substep == 0) && !(frameId % physicsParams().ccdBVHRebuildSteps);
				TICK(timeCsmpUpdatingCollisionInfoCCD);
				if (collisionParams().allowCCD)
				{
					updateCCD_CPU(rebuildCCDBVH);
				}
				TOCK_STRUCT(timeStatistics, timeCsmpUpdatingCollisionInfoCCD);

				bool rebuildDCDBVH = (substep == 0) && !(frameId % physicsParams().dcdTetMeshSceneBVHRebuildSteps);
				TICK(timeCsmpUpdatingCollisionInfoDCD);
				if (collisionParams().allowDCD)
				{
					updateDCD_CPU(rebuildDCDBVH);
				}
				TOCK_STRUCT(timeStatistics, timeCsmpUpdatingCollisionInfoDCD);

			}
			TOCK_STRUCT(timeStatistics, timeCsmpUpdatingCollisionInfoDCD);

			TICK(timeCsmpCollisionSolve);
			solveCollisionConstraintsCPU();
			TOCK_STRUCT(timeStatistics, timeCsmpCollisionSolve);

			applyPostSolvingDeformers();

		}


		// collision solve

		// run on CPU
		// syncAllToCPU();
		// std::cout << tMeshes[0]->vertPos().transpose() << std::endl;

		updateVelocities();
		++collisionDetectionCounter;

		curTime += dt;
	}
	TOCK_STRUCT(timeStatistics, timeCsmpAllSubSteps);
}

void GAIA::PBDPhysics::materialBoundarySolveGPU()
{
	bool allDone = false;
	int iParallelizationGroup = 0;
	while (!allDone)
	{
		allDone = true;
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			if (iParallelizationGroup < tMeshes[iMesh]->numAllParallelizationGroups())
			{
				//printf("---------------------------\nparallelizationGroup:\n", iParallelizationGroup);
				tMeshes[iMesh]->solveMaterialConstraintGPU_ParallelizationGroup(iParallelizationGroup, 
					physicsParams().numThreadsMaterialSolve, cudaStreams[iMesh]);
				allDone = false;

			}
		}

		++iParallelizationGroup;
	}

	//for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	//{
	//	tMeshes[iMesh]->solveMaterialGPUAggregated(physicsParams().numThreadsMaterialSolve, cudaStreams[iMesh]);
	//}

#ifdef EVAL_TET_ACCESS_COUNT
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
			cudaStreamSynchronize(cudaStreams[iMesh]);
			tMeshes[iMesh]->syncToCPU();
	}
#endif // EVAL_TET_ACCESS_COUNT

	BoundarySolve();
}

void GAIA::PBDPhysics::materialBoundarySolveGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU)
{
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->solveMaterialConstraintGPU_debugOnCPU(gpuTMeshesOnCPU[iMesh]);
	}

	BoundarySolveSolveGPU_debugOnCPU(gpuTMeshesOnCPU);
}

void GAIA::PBDPhysics::InversionSolveGPU()
{
	// solve inversions
	// verticesInversionSign has been initialized in BoundarySolve()

	for (size_t iSolve = 0; iSolve < physicsParams().inversionSolveIterations; iSolve++)
	{
		bool allDone = false;
		int iParallelizationGroup = 0;
		while (!allDone)
		{
			allDone = true;
			for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
			{
				PBDTetMeshFEM* pTM = tMeshes[iMesh].get();
				if (iParallelizationGroup < pTM->numTetParallelizationGroups())
				{
					//printf("---------------------------\nparallelizationGroup:\n", iParallelizationGroup);
					solveInversionConstraintGPU(pTM->getTetParallelizationGroupBufferGPU(iParallelizationGroup), pTM->getTetParallelizationGroupSize(iParallelizationGroup),
						physicsParams().numThreadsInversionSolve, cudaStreams[iMesh], GPUTMeshes[iMesh], physicsParams().inversionSolveConstraintMultiplier);
					allDone = false;
				}
			}

			++iParallelizationGroup;
		}

		//cpu_parallel_for(0, tMeshes.size(), [&](int iMesh) {
		//	PBDTetMeshFEM* pTM = tMeshes[iMesh].get();
		//	for (int iGroup = 0; iGroup < pTM->numTetParallelizationGroups(); iGroup++)
		//	{
		//		solveInversionConstraintGPU(pTM->getTetParallelizationGroupBufferGPU(iGroup), pTM->getTetParallelizationGroupSize(iGroup),
		//			physicsParams().numThreadsInversionSolve, cudaStreams[iMesh], GPUTMeshes[iMesh], physicsParams().inversionSolveConstraintMultiplier);
		//	}
		//	
		//});
	}
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		evaluateInversionGPU(physicsParams().numThreadsInversionSolve, cudaStreams[iMesh], tMeshes[iMesh]->numTets(), GPUTMeshes[iMesh]);
	}
}

void GAIA::PBDPhysics::BoundarySolve()
{
	auto& worldBounds = physicsParams().worldBounds;

	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		//tMeshes[iMesh]->solveMaterialConstraintGPU_ParallelizationGroup(iParallelizationGroup, physicsParams().numThreadsMaterialSolve, cudaStreams[iMesh]);
		if (physicsParams().usePlaneGround)
		{
			solveBoxBoundaryConstraint(cudaStreams[iMesh], tMeshes[iMesh]->numVertices(), physicsParams().numThreadsBoundarySolve, GPUTMeshes[iMesh], worldBounds(0, 0), worldBounds(0, 1),
				worldBounds(1, 0), worldBounds(1, 1), worldBounds(2, 0), worldBounds(2, 1), physicsParams().boundaryFrictionStatic, physicsParams().boundaryFrictionDynamic);
		}
	}
}

void GAIA::PBDPhysics::BoundarySolveSolveGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU)
{
	auto& worldBounds = physicsParams().worldBounds;
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		//tMeshes[iMesh]->solveMaterialConstraintGPU_ParallelizationGroup(iParallelizationGroup, physicsParams().numThreadsMaterialSolve, cudaStreams[iMesh]);
		if (physicsParams().usePlaneGround)
		{
			solveBoxBoundaryConstraintOnCPU(gpuTMeshesOnCPU[iMesh], worldBounds(0, 0), worldBounds(0, 1),
				worldBounds(1, 0), worldBounds(1, 1), worldBounds(2, 0), worldBounds(2, 1), physicsParams().boundaryFrictionStatic, physicsParams().boundaryFrictionDynamic);
		}
	}
}

void GAIA::PBDPhysics::collisionSolve()
{
	depthBasedCollisionSolve();

	TICK(timeCsmpCollisionSolve);
	solveCollisionConstraintsCPU();
	TOCK_STRUCT(timeStatistics, timeCsmpCollisionSolve);
}

void GAIA::PBDPhysics::depthBasedCollisionSolve()
{
	if ((collisionDetectionCounter >= physicsParams().collisionDetectionSubSteps || collisionDetectionCounter == 0)
		&& ((!physicsParams().doCollDetectionOnlyForFirstIteration) || (iIter == 0)))
	{
		if (collisionDetectionCounter != 0)
		{
			collisionDetectionCounter = 0;
		}

		positiveCollisionDetectionResults.clear();

		bool rebuildCCDBVH = (substep == 0) && !(frameId % physicsParams().ccdBVHRebuildSteps);
		TICK(timeCsmpUpdatingCollisionInfoCCD);
		if (collisionParams().allowCCD)
		{
			updateCCD_CPU(rebuildCCDBVH);
		}
		TOCK_STRUCT(timeStatistics, timeCsmpUpdatingCollisionInfoCCD);

		bool rebuildDCDBVH = (substep == 0) && !(frameId % physicsParams().dcdTetMeshSceneBVHRebuildSteps);
		TICK(timeCsmpUpdatingCollisionInfoDCD);
		if (collisionParams().allowDCD)
		{
			updateDCD_CPU(rebuildDCDBVH);
		}
		TOCK_STRUCT(timeStatistics, timeCsmpUpdatingCollisionInfoDCD);

	}
}

void GAIA::PBDPhysics::volumeBasedCollisionSolve()
{
	TICK(timeCsmpUpdatingBVHVolume);

	TOCK_STRUCT(timeStatistics, timeCsmpUpdatingBVHVolume);

}

void GAIA::PBDPhysics::initialCollisionsDetection()
{
	pDCD->updateBVH(RTC_BUILD_QUALITY_REFIT, RTC_BUILD_QUALITY_REFIT, false);

	// if CCD is not allowed, we do not need to do the initial collision detection
	// all of the vertices should be processed by DCD, thus we put all of them in initialPenetratedVertIds
	if (!collisionParams().allowCCD)
	{
		for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			PBDTetMeshFEM* pTM = tMeshes[iMesh].get();

			//initialPenetratedVertIds.resize(pTM->numSurfaceVerts());
			//for (size_t iV = 0; iV < pTM->numSurfaceVerts(); ++iV) {
			//	initialPenetratedVertIds[iV] = 
			//}
		}

		return;
	}

	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		initialPenetratedVertIds[iMesh].clear();
		initialPenetrationFreeVertIds[iMesh].clear();

		PBDTetMeshFEM* pTM = tMeshes[iMesh].get();
		std::vector<bool> vertsHasPenetrated(pTM->surfaceVIds().size());
		std::vector<CollisionDetectionResult>& dcdResults = dcdResultsAll[iMesh];

		cpu_parallel_for(0, pTM->surfaceVIds().size(), [&](int iSurfaceV)
			{
				int32_t surfaceVId = pTM->surfaceVIds()(iSurfaceV);
				pDCD->vertexCollisionDetection(surfaceVId, iMesh, &dcdResults[iSurfaceV]);

				if (dcdResults[iSurfaceV].numIntersections())
				{
					vertsHasPenetrated[iSurfaceV] = true;
				}
				else
				{
					vertsHasPenetrated[iSurfaceV] = false;
				}
			});

		for (size_t iV = 0; iV < vertsHasPenetrated.size(); ++iV) {
			if (vertsHasPenetrated[iV])
			{
				initialPenetratedVertIds[iMesh].push_back(iV);

			}
			else {
				initialPenetrationFreeVertIds[iMesh].push_back(iV);

			}
		}
	}
}

void GAIA::PBDPhysics::updateDCD_CPU(bool rebuild)
{
	//positiveCollisionDetectionResults.clear();
	TICK(timeCsmpUpdatingBVHDCD);
	int numActiveCCD = positiveCollisionDetectionResults.size();

	RTCBuildQuality surfaceSceneQuality 
		= !(frameId % physicsParams().dcdSurfaceSceneBVHRebuildSteps) && (substep == 0) ? RTC_BUILD_QUALITY_LOW : RTC_BUILD_QUALITY_REFIT;
	if (rebuild)
	{
		pDCD->updateBVH(RTC_BUILD_QUALITY_LOW, surfaceSceneQuality, true);
	}
	else
	{
		pDCD->updateBVH(RTC_BUILD_QUALITY_REFIT, surfaceSceneQuality, true);
	}
	TOCK_STRUCT(timeStatistics, timeCsmpUpdatingBVHDCD);

	if (!collisionParams().allowCCD)
	{
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			TetMeshFEM* pTM = tMeshes[iMesh].get();
			std::vector<CollisionDetectionResult>& dcdResults = dcdResultsAll[iMesh];

			TICK(timeCsmpColDetectDCD);
			cpu_parallel_for(0, pTM->surfaceVIds().size(), [&](int iSurfaceV)
				{
					int32_t surfaceVId = pTM->surfaceVIds()(iSurfaceV);
					if (pTM->DCDEnabled(surfaceVId))
					{
						pDCD->vertexCollisionDetection(surfaceVId, iMesh, &dcdResults[iSurfaceV]);
					}
				});
			TOCK_STRUCT(timeStatistics, timeCsmpColDetectDCD);

			// recorded the surface vert ids that are tested positive for collision
			std::vector<int> dcdResultsPositiveVId;
			for (int iSurfaceV = 0; iSurfaceV < pTM->surfaceVIds().size(); iSurfaceV++)
			{
				if (dcdResults[iSurfaceV].numIntersections())
				{
					dcdResultsPositiveVId.push_back(iSurfaceV);
				}
			}

			TICK(timeCsmpShortestPathSearchDCD);
			cpu_parallel_for(0, dcdResultsPositiveVId.size(), [&](int iIntersection)
				{
					ClosestPointQueryResult closestPtResult;
					pDCD->closestPointQuery(&dcdResults[dcdResultsPositiveVId[iIntersection]], &closestPtResult);
				});
			TOCK_STRUCT(timeStatistics, timeCsmpShortestPathSearchDCD);

			for (int iIntersection = 0; iIntersection < dcdResultsPositiveVId.size(); iIntersection++)
			{
				positiveCollisionDetectionResults.push_back(dcdResults[dcdResultsPositiveVId[iIntersection]]);
			}
		}
	}
	else
	{
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			if (!initialPenetratedVertIds[iMesh].size())
			{
				continue;
			}
			TetMeshFEM* pTM = tMeshes[iMesh].get();
			std::vector<CollisionDetectionResult> dcdResults(initialPenetratedVertIds[iMesh].size());

			TICK(timeCsmpColDetectDCD);
			cpu_parallel_for(0, dcdResults.size(), [&](int iSurfaceV)
				{
					int32_t surfaceVId = pTM->surfaceVIds()(initialPenetratedVertIds[iMesh][iSurfaceV]);
					if (pTM->DCDEnabled(surfaceVId))
					{
						pDCD->vertexCollisionDetection(surfaceVId, iMesh, &dcdResults[iSurfaceV]);
					}
				});
			TOCK_STRUCT(timeStatistics, timeCsmpColDetectDCD);

			// recorded the surface vert ids that are tested positive for collision
			std::vector<int> dcdResultsPositiveVId;
			for (int iSurfaceV = 0; iSurfaceV < dcdResults.size(); iSurfaceV++)
			{
				if (dcdResults[iSurfaceV].numIntersections())
				{
					dcdResultsPositiveVId.push_back(iSurfaceV);
				}
			}

			TICK(timeCsmpShortestPathSearchDCD);
			cpu_parallel_for(0, dcdResultsPositiveVId.size(), [&](int iIntersection)
				{
					ClosestPointQueryResult closestPtResult;
					pDCD->closestPointQuery(&dcdResults[dcdResultsPositiveVId[iIntersection]], &closestPtResult);
				});
			TOCK_STRUCT(timeStatistics, timeCsmpShortestPathSearchDCD);

			for (int iIntersection = 0; iIntersection < dcdResultsPositiveVId.size(); iIntersection++)
			{
				positiveCollisionDetectionResults.push_back(dcdResults[dcdResultsPositiveVId[iIntersection]]);
			}
		}
	}
	int numActiveDCD = positiveCollisionDetectionResults.size() - numActiveCCD;


#ifdef DO_COLLISION_STATISTICS
	collisionStatistics.numOfCollisionsDCDs.push_back(numActiveDCD);
	collisionStatistics.numOfCollisionsCCDs.push_back(numActiveCCD);
	collisionStatistics.numOfBVHQuerysEachStep.emplace_back();
	collisionStatistics.numOfTetTraversal.emplace_back();
	collisionStatistics.numTetsTraversed.emplace_back();
	for (int iCol = numActiveCCD; iCol < positiveCollisionDetectionResults.size(); iCol++)
	{
		collisionStatistics.numOfBVHQuerysEachStep.back().push_back(positiveCollisionDetectionResults[iCol].numberOfBVHQuery);
		collisionStatistics.numOfTetTraversal.back().push_back(positiveCollisionDetectionResults[iCol].numberOfTetTraversal);
		collisionStatistics.numTetsTraversed.back().push_back(positiveCollisionDetectionResults[iCol].numberOfTetsTraversed);
	}
#endif // DO_COLLISION_STATISTICS

	// std::cout << "Number of collisions detected by DCD:" << numActiveDCD << std::endl;

}

void GAIA::PBDPhysics::updateCCD_CPU(bool rebuild)
{
	TICK(timeCsmpUpdatingBVHCCD);
	if (rebuild)
	{
		pCCD->updateBVH(RTC_BUILD_QUALITY_LOW);
	}
	else
	{
		pCCD->updateBVH(RTC_BUILD_QUALITY_REFIT);
	}
	TOCK_STRUCT(timeStatistics, timeCsmpUpdatingBVHCCD);

	int numActiveCCD = positiveCollisionDetectionResults.size();

	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		TetMeshFEM* pTM = tMeshes[iMesh].get();

		TICK(timeCsmpColDetectCCD);
		cpu_parallel_for(0, initialPenetrationFreeVertIds[iMesh].size(), [&](int iSurfaceV)
			{
				int32_t surfaceVIdSurfaceMesh = initialPenetrationFreeVertIds[iMesh][iSurfaceV];
				int32_t surfaceVId = pTM->surfaceVIds()(surfaceVIdSurfaceMesh);
				if (pTM->CCDEnabled(surfaceVId))
				{
					pCCD->vertexContinuousCollisionDetection(surfaceVId, iMesh, &ccdResultsAll[iMesh][surfaceVIdSurfaceMesh]);
				}
			});

		for (int iSurfaceV = 0; iSurfaceV < initialPenetrationFreeVertIds[iMesh].size(); iSurfaceV++)
		{
			int32_t surfaceVIdSurfaceMesh = initialPenetrationFreeVertIds[iMesh][iSurfaceV];
			if (ccdResultsAll[iMesh][surfaceVIdSurfaceMesh].numIntersections())
			{
				positiveCollisionDetectionResults.push_back(ccdResultsAll[iMesh][surfaceVIdSurfaceMesh]);
			}
		}
		TOCK_STRUCT(timeStatistics, timeCsmpColDetectCCD);
	}
	numActiveCCD = positiveCollisionDetectionResults.size() - numActiveCCD;

	// std::cout << "Number of collisions detected by CCD:" << numActiveCCD << std::endl;

}

void GAIA::PBDPhysics::changeNumSubSteps(int newNumSubSteps)
{
	dt = physicsAllParams.physicsParams.timeStep / physicsAllParams.physicsParams.numSubsteps;

	// make parameters step invariant
	for (size_t iMesh = 0; iMesh < objectParamsList.objectParams.size(); iMesh++)
	{
		objectParamsList.getObjectParam(iMesh).exponentialVelDamping = std::pow(objectParamsList.getObjectParam(iMesh).exponentialVelDamping, dt);
		objectParamsList.getObjectParam(iMesh).constantVelDamping = objectParamsList.getObjectParam(iMesh).constantVelDamping * dt;
	}
	for (int iMesh = 0; iMesh < tMeshes.size(); ++iMesh)
	{
		GAIA::PBDTetMeshFEM* pTM = tMeshes[iMesh].get();

		GAIA::PBDTetMeshNeoHookean* pTMNeoHookean;
		GAIA::PBDTetMeshMassSpring* pTMMassSpring;
		switch (pTM->pObjectParams->materialType)
		{
		case GAIA::NeoHookean:
			pTMNeoHookean = (GAIA::PBDTetMeshNeoHookean*)pTM;
			pTMNeoHookean->dt = dt;
			pTMNeoHookean->tetMeshGPU.dt = dt;
			pTMNeoHookean->tetMeshGPU_forCPU.dt = dt;
		default:
			pTMMassSpring = (GAIA::PBDTetMeshMassSpring*)pTM;
			pTMMassSpring->dt = dt;
			pTMMassSpring->tetMeshGPU.dt = dt;
			pTMMassSpring->tetMeshGPU_forCPU.dt = dt;
			break;
		}
	}

	updateGPUMeshes();
}

void GAIA::PBDPhysics::updateVelocities()
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		tMeshes[iMesh]->mVelocity = (tMeshes[iMesh]->mVertPos - tMeshes[iMesh]->mVertPrevPos) * (1.0f / dt);

		VecDynamic velMags = tMeshes[iMesh]->mVelocity.colwise().norm();
		FloatingType maxVelocityMagnitude = getObjectParam(iMesh).maxVelocityMagnitude;

		//for (size_t iVert = 0; iVert < tMeshes[iMesh]->numVertices(); iVert++)
		cpu_parallel_for(0, tMeshes[iMesh]->numVertices(), [&](int iVert)
			{
				FloatingType vMag = velMags(iVert);

				if (getObjectParam(iMesh).hasNoGravZone &&
					maxVelocityMagnitude > 0 &&
					tMeshes[iMesh]->mVertPos(GRAVITY_AXIS, iVert) > getObjectParam(iMesh).noGravZoneThreshold)
					// no vel damping in no gravity zone
					// apply the maximum velocity constraint
				{
					if (vMag > maxVelocityMagnitude) {
						tMeshes[iMesh]->mVelocity.col(iVert) *= (maxVelocityMagnitude / vMag);
					}
				}
				else if (vMag > 1e-6f) {
					FloatingType vMagNew;
					vMagNew = vMag * getObjectParam(iMesh).exponentialVelDamping - getObjectParam(iMesh).constantVelDamping;
					vMagNew = vMagNew > 1e-6f ? vMagNew : 0.f;
					tMeshes[iMesh]->mVelocity.col(iVert) *= vMagNew / vMag;
				}
				else
				{
					tMeshes[iMesh]->mVelocity.col(iVert) *= 0.f;
				}
		});

	}
}

void GAIA::PBDPhysics::updateWorldBox()
{
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		bool changed = false;
		TetMeshFEM* pTM = tMeshes[iMesh].get();
		for (int iDim = 0; iDim < 3; iDim++)
		{
			FloatingType dimMin = pTM->vertices().row(iDim).minCoeff();
			FloatingType dimMax = pTM->vertices().row(iDim).maxCoeff();
			if (physicsParams().worldBounds(iDim, 0) > dimMin)
			{
				physicsParams().worldBounds(iDim, 0) = dimMin;
				changed = true;
			}

			if (physicsParams().worldBounds(iDim, 1) < dimMax)
			{
				physicsParams().worldBounds(iDim, 1) = dimMax;
				changed = true;
			}
		}

		std::cout << "World box updates to: \n" << physicsParams().worldBounds.transpose() << std::endl;
	}
}


void GAIA::PBDPhysics::solveCollisionConstraintsCPU()
{
	int numCollisions = positiveCollisionDetectionResults.size();

	for (int iPenetratedV = 0; iPenetratedV < numCollisions; iPenetratedV++)
	{
		int curVertID = positiveCollisionDetectionResults[iPenetratedV].idVQuery;
		int curMeshID = positiveCollisionDetectionResults[iPenetratedV].idTMQuery;

		CollisionDetectionResult & colResult = positiveCollisionDetectionResults[iPenetratedV];

		for (int iIntersection = 0; iIntersection < colResult.numIntersections(); iIntersection++)
		{
			if (colResult.collidingPts[iIntersection].shortestPathFound) {
				solveCollisionConstraintsOneVertCPU(curVertID, curMeshID, colResult, iIntersection);
			}
		}
	}
}


void GAIA::PBDPhysics::solveCollisionConstraintsOneVertCPU(int32_t curVertID, int32_t curMeshID,
	CollisionDetectionResult& colResult, int32_t iIntersection)
{
	int collidedMeshID = colResult.collidingPts[iIntersection].intersectedMeshId;
	int closestFaceId = colResult.collidingPts[iIntersection].closestSurfaceFaceId;

	PBDTetMeshFEM* pCurTM = tMeshes[curMeshID].get();
	PBDTetMeshFEM* pIntersectedTM = tMeshes[collidedMeshID].get();


	if (collidedMeshID < 0 && closestFaceId < 0) {
		std::cout << "CLOSEST TRI-ID is NEGATIVE!!!" << std::endl;

		return;
	}

	ClosestPointOnTriangleType pointType;

	Vec3 curVertPos = pCurTM->mVertPos.col(colResult.idVQuery);
	Vec3 closestSurfacePoint, barycentrics;
	Vec3I closestFaceVIds = pIntersectedTM->surfaceFacesTetMeshVIds().col(closestFaceId);

	if (
		//(iIter == 0 || !physicsParams().doCollDetectionOnlyForFirstIteration) 
		//&& !colResult.fromCCD
		// collision solution may influence the shape of the mesh thus the closest point information maybe inaccurate even for the first iteration
		collisionParams().restPoseCloestPoint
		)
		// the closest point is still valid
	{
		pointType = colResult.collidingPts[iIntersection].closestPointType;
		closestSurfacePoint = colResult.collidingPts[iIntersection].closestSurfacePt;
		barycentrics = colResult.collidingPts[iIntersection].closestSurfacePtBarycentrics;
	}
	else
		// recompute the closest point
	{
		embree::Vec3fa a = loadVertexPos(pIntersectedTM, pIntersectedTM->surfaceFacesTetMeshVIds()(0, closestFaceId)),
			b = loadVertexPos(pIntersectedTM, pIntersectedTM->surfaceFacesTetMeshVIds()(1, closestFaceId)),
			c = loadVertexPos(pIntersectedTM, pIntersectedTM->surfaceFacesTetMeshVIds()(2, closestFaceId)),
			p = loadVertexPos(pCurTM, colResult.idVQuery);

		embree::Vec3fa bary;
		embree::Vec3fa closestSurfacePoint_ = closestPointTriangle(p, a, b, c, bary, pointType);

		closestSurfacePoint << closestSurfacePoint_.x, closestSurfacePoint_.y, closestSurfacePoint_.z;
		barycentrics << bary.x, bary.y, bary.z;
	}
	
	Vec3 normal;
	int32_t surfaceVIdTMeshIndex = -1;
	int32_t surfaceVIdSurfaceIndex = -1;
	switch (pointType)
	{
	case GAIA::ClosestPointOnTriangleType::AtA:
		surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(0, closestFaceId);
		surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(0, closestFaceId);
		pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);

		break;
	case GAIA::ClosestPointOnTriangleType::AtB:
		surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(1, closestFaceId);
		surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(1, closestFaceId);
		pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::AtC:
		surfaceVIdTMeshIndex = pIntersectedTM->surfaceFacesTetMeshVIds()(2, closestFaceId);
		surfaceVIdSurfaceIndex = pIntersectedTM->surfaceFacesSurfaceMeshVIds()(2, closestFaceId);
		pIntersectedTM->computeVertexNormal(surfaceVIdSurfaceIndex, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::AtAB:
		pIntersectedTM->computeEdgeNormal(closestFaceId, 0, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::AtBC:
		pIntersectedTM->computeEdgeNormal(closestFaceId, 1, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::AtAC:
		pIntersectedTM->computeEdgeNormal(closestFaceId, 2, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::AtInterior:
		pIntersectedTM->computeFaceNormal(closestFaceId, normal);
		break;
	case GAIA::ClosestPointOnTriangleType::NotFound:
		return;
		break;
	default:
		return;
		break;
	}

	Vec3 r = curVertPos - closestSurfacePoint;

	// when models are pressed to corners and has degenerative faces
	FloatingType nNorm = normal.norm();
	FloatingType penetrationDepth = r.norm();
	//FloatingType penetrationDepth = -((r * n) / n.norm());

	if (isnan(nNorm) || nNorm == 0.f)
	{
		return;
	}
	else if (r.norm() == 0.f)
	{
		return;
	}

	FloatingType C = r.dot( normal);
	// collision constraint has been met
	if (C > 0.f)
	{
		C = 0.f;
	}
	Vec3 dC0 = normal * (-barycentrics(0));
	Vec3 dC1 = normal * (-barycentrics(1));
	Vec3 dC2 = normal * (-barycentrics(2));

	Vec12 grads;
	grads << dC0[0]
		, dC0[1]
		, dC0[2]
		, dC1[0]
		, dC1[1]
		, dC1[2]
		, dC2[0]
		, dC2[1]
		, dC2[2]
		, normal[0]
		, normal[1]
		, normal[2];

	applyCollisionConstraintsCPU(closestFaceVIds, curVertID, curMeshID, collidedMeshID, C, penetrationDepth, dt, barycentrics, grads);
}

void GAIA::PBDPhysics::applyCollisionConstraintsCPU(const Vec3I& closestFaceVIds, int32_t curVertId, int32_t curMeshID, int32_t collidedMeshID, FloatingType C,
	FloatingType penetrationDepth, FloatingType dt, Vec3& barycentric, const Vec12& gradients)
{
	C = C * physicsParams().collisionConstraintMultiplier;

	if (C == 0.0) {
		return;
	}

	PBDTetMeshFEM* pCurTM = tMeshes[curMeshID].get();
	PBDTetMeshFEM* pIntersectedTM = tMeshes[collidedMeshID].get();

	FloatingType w = 0.0;
	FloatingType damp = 0.0;

	for (int i = 0; i < 3; i++) {
		int vId = closestFaceVIds[i];
		
		w += gradients.block<3, 1>(i * 3, 0).squaredNorm() * pIntersectedTM->vertexInvMass(vId); // ((*softBodies)[mID]->invMass[id]);
	}

	w += gradients.block<3, 1>(9, 0).squaredNorm() * pCurTM->vertexInvMass(curVertId);

	if (w == 0.0) {
		return;
	}

	FloatingType dlambda = -C / (w);

	std::array<Vec3, 4>  velocities;
	for (int i = 0; i < 3; i++) {
		int vId = closestFaceVIds[i];
		Vec3 correction = dlambda * pIntersectedTM->vertexInvMass(vId) * gradients.block<3, 1>(i * 3, 0);
		pIntersectedTM->vertex(vId) += correction;
		velocities[i] = pIntersectedTM->vertex(vId) - pIntersectedTM->vertexPrevPos(vId);
	}

	Vec3 correction = dlambda * pCurTM->vertexInvMass(curVertId) * gradients.block<3, 1>(9, 0);
	pCurTM->vertex(curVertId) += correction;
	velocities[3] = pCurTM->vertex(curVertId) - pCurTM->vertexPrevPos(curVertId);

	Vec3 g3N = gradients.block<3, 1>(9, 0);
	g3N /= g3N.norm();

	Vec3 closestPVelocity = velocities[0] * barycentric[0] + velocities[1] * barycentric[1] + velocities[2] * barycentric[2];
	Vec3 perpVelocity = closestPVelocity - g3N * closestPVelocity.dot(g3N);

	Vec3 collidedPointPerpVelocity = velocities[3] - g3N * velocities[3].dot(g3N);

	Vec3 diff = collidedPointPerpVelocity - perpVelocity;

	FloatingType dNorm = diff.norm();
	FloatingType dynamicNorm;
	FloatingType frictionRatioDynamic = objectParamsList.getObjectParam(curMeshID).friction_dynamic;
	FloatingType frictionRatioStatic = objectParamsList.getObjectParam(curMeshID).friction_static;

	if (penetrationDepth * frictionRatioDynamic / dNorm <= 1)
	{
		dynamicNorm = penetrationDepth * frictionRatioDynamic;
	}
	else
	{
		dynamicNorm = dNorm;
	}

	dynamicNorm = penetrationDepth * frictionRatioDynamic;


	FloatingType staticNorm = penetrationDepth * frictionRatioStatic;

	if (dNorm > staticNorm) {
		diff = (diff / dNorm) * dynamicNorm;
	}

	Vec3 diffV = diff * pCurTM->vertexInvMass(curVertId) / (pCurTM->vertexInvMass(curVertId)
		+ pIntersectedTM->vertexInvMass(closestFaceVIds[0]) 
		+ pIntersectedTM->vertexInvMass(closestFaceVIds[1]) 
		+ pIntersectedTM->vertexInvMass(closestFaceVIds[2]));
	pCurTM->vertex(curVertId) -= diff ;

	Vec3 diffTriangle = diff - diffV;

	//FloatingType d = 0.f;
	//for (int i = 0; i < 3; i++) {
	//	int vId = closestFaceVIds[i];
	//	dlambda* pIntersectedTM->vertexInvMass(vId);
	//}

	//FloatingType d = pIntersectedTM->vertexInvMass(closestFaceVIds[0]) * pow(barycentric[0],2)
	//	+ pIntersectedTM->vertexInvMass(closestFaceVIds[1]) * pow(barycentric[1], 2)
	//	+ pIntersectedTM->vertexInvMass(closestFaceVIds[2]) * pow(barycentric[2], 2);

	//pIntersectedTM->vertex(closestFaceVIds[0]) += diffTriangle * pIntersectedTM->vertexInvMass(closestFaceVIds[0]) * barycentric[0] / d;
	//pIntersectedTM->vertex(closestFaceVIds[1]) += diffTriangle * pIntersectedTM->vertexInvMass(closestFaceVIds[1]) * barycentric[1] / d;
	//pIntersectedTM->vertex(closestFaceVIds[2]) += diffTriangle * pIntersectedTM->vertexInvMass(closestFaceVIds[2]) * barycentric[2] / d;

	pIntersectedTM->vertex(closestFaceVIds[0]) += diffTriangle ;
	pIntersectedTM->vertex(closestFaceVIds[1]) += diffTriangle ;
	pIntersectedTM->vertex(closestFaceVIds[2]) += diffTriangle ;

	//for (int i = 0; i < 3; i++) {
	//	int vId = closestFaceVIds[i];
	//	Vec3 diffPerVert = dlambda * pIntersectedTM->vertexInvMass(vId) * gradients.block<3, 1>(i * 3, 0);
	//	pIntersectedTM->vertex(vId) += correction;
	//	velocities[i] = pIntersectedTM->vertex(vId) - pIntersectedTM->vertexPrevPos(vId);
	//}

	//Vec3 correction = dlambda * pCurTM->vertexInvMass(curVertId) * gradients.block<3, 1>(9, 0);
	//pCurTM->vertex(curVertId) += correction;
}

void GAIA::PBDPhysics::writeOutputs(std::string outFolder, int frameId)
{
	std::ostringstream aSs;
	aSs << std::setfill('0') << std::setw(8) << frameId;
	std::string outNumber = aSs.str();
	std::string outFile = outFolder + "/A" + outNumber + "." + physicsAllParams.physicsParams.outputExt;

	std::string tetMeshOutOutPath = outFolder + "/TetMesh";
	std::string tetMeshOutStatistics = outFolder + "/Statistics";
	std::vector<TetMeshFEM::SharedPtr> basetetMeshes;
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		basetetMeshes.push_back(std::static_pointer_cast<TetMeshFEM>(tMeshes[iMesh]));
	}
	if (physicsAllParams.physicsParams.outputExt == "bin")
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

		if (!(frameId % physicsAllParams.physicsParams.binaryModeVisualizationSteps))
		{
			std::string outFileVis = outFolder + "/A" + outNumber + ".ply";

			writeAllToPLY(outFileVis.c_str(), basetetMeshes, physicsAllParams.physicsParams.saveAllModelsTogether, physicsAllParams.physicsParams.saveAllModelsTogether);
		}
	}
	else if (physicsAllParams.physicsParams.outputExt == "ply")
	{
		
		writeAllToPLY(outFile.c_str(), basetetMeshes, physicsAllParams.physicsParams.saveAllModelsTogether);
	}
	else if (physicsAllParams.physicsParams.outputExt == "obj") {
		// writeAllToObj(outFile.c_str(), getSoftBodies(), physicsAllParams.physicsParams.saveAllModelsTogether);
	}
	else
	{
		std::cout << "[Error] Unsupported output EXT name: " << physicsAllParams.physicsParams.outputExt << std::endl;
		return;
	}

	if (physicsAllParams.physicsParams.outputT) {
		std::string outFileT = tetMeshOutOutPath + "/A" + outNumber + ".t";
		//writeAllTetMeshToT(outFileT.c_str(), getSoftBodies());
	}

	if (physicsAllParams.physicsParams.outputVTK) {
		std::string outFileVtk = tetMeshOutOutPath + "/A" + outNumber + ".vtk";
		// writeAllTetMeshToVtk(outFileVtk.c_str(), getSoftBodies());
	}

	if (physicsAllParams.physicsParams.outputStatistics)
	{
		//pSoftbodyManager->statistics.writeToJsonFile(outFolder + "/Statistics.json");
		timeStatistics.writeToJsonFile(tetMeshOutStatistics + "/A" + outNumber + ".json");
		if (frameId % physicsAllParams.physicsParams.outputRecoveryStateStep == 0) {
#ifdef DO_COLLISION_STATISTICS
			collisionStatistics.writeToJsonFile(tetMeshOutStatistics + "/CollisionStatisticsAll.json");
#endif // DO_COLLISION_STATISTICS

		}
	}

	if (physicsAllParams.physicsParams.outputRecoveryState && !(frameId % physicsAllParams.physicsParams.outputRecoveryStateStep) && (frameId))
	{
		std::string stateOutOutPath = outFolder + "/RecoveryStates";
		PBDPhysicsState state;
		state.fromPhysics(*this, frameId);
		std::string outFileState = stateOutOutPath + "/A" + outNumber + ".json";
		state.writeToJsonFile(outFileState);
	}

}


bool GAIA::PBDPhysics::initializeFromState(std::string& stateFile)
{
	// PBDPhysicsState state;
	std::cout << "----------------------------------------------------\n"
		<< "Recovering state from:" << stateFile << "\n"
		<< "----------------------------------------------------\n";
	// state.loadFromJsonFile(stateFile);
	// state.initializeSoftBodyManager(*pSoftbodyManager, pSoftbodyManager->frameId);

	return false;
}

void GAIA::PBDPhysics::setUpOutputFolders(std::string outFolder)
{
	MF::IO::createFolder(outFolder);

	std::string tetMeshOutOutPath = outFolder + "/TetMesh";
	if (physicsAllParams.physicsParams.outputVTK || physicsAllParams.physicsParams.outputT)
	{
		MF::IO::createFolder(tetMeshOutOutPath);
	}

	if (physicsAllParams.physicsParams.saveSimulationParameters)
	{
		std::string paramsOutOutPath = outFolder + "/SimulationParameters_" + MF::Time::getCurTimeString();
		MF::IO::createFolder(paramsOutOutPath);
		saveExperimentParameters(paramsOutOutPath);
	}

	if (physicsAllParams.physicsParams.outputStatistics) {
		std::string statisticsOutOutPath = outFolder + "/Statistics";
		MF::IO::createFolder(statisticsOutOutPath);
	}

	if (physicsAllParams.physicsParams.outputRecoveryState) {
		std::string stateOutOutPath = outFolder + "/RecoveryStates";
		MF::IO::createFolder(stateOutOutPath);
	}

}

void GAIA::PBDPhysics::applyInitialGuesses()
{
	cpu_parallel_for(0, tMeshes.size(), [&] (int iMesh) {
		tMeshes[iMesh]->applyInitialGuess();
		});
}

void GAIA::PBDPhysics::evaluateObjectActiveness()
{
}

void GAIA::PBDPhysics::applyDeformers()
{
	// apply deformers
	for (size_t iDeformer = 0; iDeformer < deformers.size(); iDeformer++)
	{
		(*deformers[iDeformer])(this, curTime, frameId, substep, iIter, dt);
	}
}

void GAIA::PBDPhysics::applyConstraints()
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		PBDTetMeshFEM* pTM = tMeshes[iMesh].get();
		for (size_t iFixPoint = 0; iFixPoint < pTM->pObjectParams->fixedPoints.size(); iFixPoint++)
		{
			int vertexId = pTM->pObjectParams->fixedPoints[iFixPoint];
			pTM->vertex(vertexId) = pTM->mVertPrevPos.col(vertexId);
		}
	}
}

void GAIA::PBDPhysics::applyPostSolvingDeformers()
{
	// apply deformers
	for (size_t iDeformer = 0; iDeformer < postSolvingDeformers.size(); iDeformer++)
	{
		(*postSolvingDeformers[iDeformer])(this, curTime, frameId, substep, iIter, dt);
	}
}

void GAIA::PBDPhysics::solveInvertedTets()
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		PBDTetMeshFEM* pTM = tMeshes[iMesh].get();

		cpu_parallel_for(0, pTM->numTets(), [&](int iTet) {
			FloatingType v = CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data());
			if (v <= 0)
			{
				pTM->tetsInvertedSign[iTet] = true;
			}
			else
			{
				pTM->tetsInvertedSign[iTet] = false;
			}
			});

		//int numInverted = 0;
		//for (int iTet = 0; iTet < pTM->numTets(); iTet++)
		//{
		//	if (CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data()) < 0) {
		//		++numInverted;
		//	}
		//}
		//std::cout << "num inverted before solving: " << numInverted << std::endl; 
		if (physicsParams().solveInvertedTets)
		{
			// -- v2
			const int inversionSolveIters = 2;
			for (int i_ = 0; i_ < inversionSolveIters; i_++)
			{
				//std::vector<int> checked(pTM->numTets(), false);
				for (int iTet = 0; iTet < pTM->numTets(); iTet++)
				{
					if (pTM->tetsInvertedSign[iTet])
					{
						solveInversionConstraint(pTM, iTet);
						if (CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data()) > 0)
						{
							pTM->tetsInvertedSign[iTet] = false;
						}

#ifdef TET_TET_ADJACENT_LIST
						//for (int iNeiTet = 0; iNeiTet < pTM->tetAllNeighborTets[iTet].size(); iNeiTet++)
						cpu_parallel_for(0, pTM->tetAllNeighborTets[iTet].size(), [&](int iNeiTet)
							{
								int neiTetId = pTM->tetAllNeighborTets[iTet][iNeiTet];
								if (!pTM->tetsInvertedSign[neiTetId])
									if ( CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(neiTetId).data()) <= 0)
									{
										pTM->tetsInvertedSign[neiTetId] = true;
									}
								//checked[neiTetId] = true;
							});
#endif // TET_TET_ADJACENT_LIST

					}
				}
			}

			//numInverted = 0;
			//for (int iTet = 0; iTet < pTM->numTets(); iTet++)
			//{
			//	if (CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data()) < 0) {
			//		++numInverted;
			//	}
			//}
			//std::cout << "num inverted after solving: " << numInverted << std::endl;

			// -- v1
			//std::unordered_set<int32_t> possibleInvertedTetIds;
			//int numInverted = 0;
			//for (int iTet = 0; iTet < pTM->numTets(); iTet++)
			//{
			//	if (pTM->tetsInvertedSign[iTet])
			//	{
			//		solveInversionConstraint(pTM, iTet);
			//		possibleInvertedTetIds.insert(iTet);

			//		possibleInvertedTetIds.insert(pTM->tetAllNeighborTets[iTet].cbegin(), (pTM->tetAllNeighborTets[iTet].cend()));
			//		if ( CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data()) > 0)
			//		{
			//			pTM->tetsInvertedSign[iTet] = false;
			//		}
			//		++numInverted;
			//	}
			//}

			//int numAttempts = 0;
			//while (!possibleInvertedTetIds.empty())
			//{
			//	auto pBegin = possibleInvertedTetIds.begin();
			//	int iTet = *pBegin;
			//	bool isSurfaceTet = pTM->tetsIsSurfaceTet[iTet];

			//	if (CuMatrix::tetOrientedVolume(pTM->mVertPos.data(), pTM->tetVIds().col(iTet).data()) > 0)
			//	{
			//		possibleInvertedTetIds.erase(pBegin);
			//		pTM->tetsInvertedSign[iTet] = false;
			//	}
			//	else
			//	{
			//		if (numAttempts <= physicsParams().maxInversionSolveMultiplier * numInverted || isSurfaceTet)
			//		{
			//			++numAttempts;
			//			solveInversionConstraint(pTM, iTet);
			//			possibleInvertedTetIds.insert(pTM->tetAllNeighborTets[iTet].cbegin(), pTM->tetAllNeighborTets[iTet].cend());
			//		}
			//		else
			//			// remove non-surface tets when numAttempts > physicsParams().maxInversionSolveMultiplier
			//		{
			//			possibleInvertedTetIds.erase(pBegin);
			//		}

			//		if (numAttempts > physicsParams().maxSurfaceInversionSolveMultiplier * numInverted)
			//		{
			//			break;
			//		}
			//	}

			//}
		}

		//std::fill(pTM->verticesInvertedSign.begin(), pTM->verticesInvertedSign.end(), false);
		pTM->verticesInvertedSign.fill(false);
		//for (int iTet = 0; iTet < pTM->numTets(); iTet++)
		cpu_parallel_for(0, pTM->numTets(), [&](int iTet)
			{
				if (pTM->tetsInvertedSign[iTet])
				{
					for (size_t iV = 0; iV < 4; iV++)
					{
						pTM->verticesInvertedSign[pTM->getTVId(iTet, iV)] = true;
					}
				}
			});
	}
}

void GAIA::PBDPhysics::solveInversionConstraint(PBDTetMeshFEM* pTM, int iTet)
{
	int id0 = pTM->tetVIds()(0, iTet);
	int id1 = pTM->tetVIds()(1, iTet);
	int id2 = pTM->tetVIds()(2, iTet);
	int id3 = pTM->tetVIds()(3, iTet);

	Vec3 e1 = pTM->vertex(id1) - pTM->vertex(id0);
	Vec3 e2 = pTM->vertex(id2) - pTM->vertex(id0);
	Vec3 e3 = pTM->vertex(id3) - pTM->vertex(id0);

	double vol = (e1.dot(e2.cross(e3))) / 6.0f;

	Vec3 gradC0 = (e3 - e1).cross(e2 - e1);
	Vec3 gradC1 = (e2.cross(e3));
	Vec3 gradC2 = (e3.cross(e1));
	Vec3 gradC3 = (e1.cross(e2));

	std::array<FloatingType, 12> grads =
	{
		gradC0[0],
		gradC0[1],
		gradC0[2],
		gradC1[0],
		gradC1[1],
		gradC1[2],
		gradC2[0],
		gradC2[1],
		gradC2[2],
		gradC3[0],
		gradC3[1],
		gradC3[2]
	};

	if (vol < 0) {
		FloatingType C = 6 * (vol - (1.f / pTM->tetInvRestVolume(iTet)) * physicsParams().inversionSolveConstraintMultiplier);
		applyToInfiniteStiffness(pTM->mVertPos.data(), pTM->vertexInvMass.data(), iTet, pTM->tetVIds().col(iTet).data(), C, grads.data());
	}
}

void GAIA::PBDPhysics::evaluateConvergence(int iStep, int iIteration)
{
	if (!physicsParams().evaluateConvergence)
	{
		return;
	}
	if (timeStatistics.meritEnergy.size() < iStep + 1)
	{
		timeStatistics.meritEnergy.emplace_back();
	}
	timeStatistics.meritEnergy.back().emplace_back();
	for (size_t iObj = 0; iObj < tMeshes.size(); iObj++)
	{
		timeStatistics.meritEnergy[iStep].back().push_back(tMeshes[iObj]->evaluateMeritEnergy());
	}

}

std::vector<GAIA::DeformerPtr>& GAIA::PBDPhysics::getDeformers()
{
	return deformers;
}

std::vector<GAIA::DeformerPtr>& GAIA::PBDPhysics::getPostSolvingDeformers() {
	return postSolvingDeformers;
};


size_t GAIA::ObjectsParamList::size()
{
	return objectParams.size();
}


bool GAIA::ObjectsParamList::fromJson(nlohmann::json& objectParam)
{
	nlohmann::json modelsInfo = objectParam["Models"];

	for (nlohmann::json& model : modelsInfo) {
		std::string materialName;
		parseJsonParameters(model, "materialName", materialName);

		if (materialName == "NeoHookean")
		{
			objectParams.push_back(std::make_shared<ObjectParamsPBDNeoHookean>());
			objectParams.back()->fromJson(model);
			objectParams.back()->materialName = materialName;
			// tMeshes.push_back(std::make_shared<PBDTetMeshNeoHookean>());
		}
		else if (materialName == "MassSpring") {
			// wait to be fullfilled
			objectParams.push_back(std::make_shared<ObjectParamsPBDMassSpring>());
			objectParams.back()->fromJson(model);
			objectParams.back()->materialName = materialName;
		}
		else
		{
			std::cout << "Warning!!! Material name: " << materialName << " not recognized! Skipping this material!\n";
			/*std::cout << "Using default material NeoHookean instead!\n";
			objectParams.push_back(std::make_shared<ObjectParamsPBDNeoHookean>());*/
		}
	}

	return true;
}

bool GAIA::ObjectsParamList::toJson(nlohmann::json& objectParam)
{
	for (size_t iObj = 0; iObj < objectParams.size(); iObj++)
	{
		nlohmann::json obj;
		objectParams[iObj]->toJson(obj);
		objectParam["Models"].push_back(obj);
	}
	return true;
}

bool GAIA::PBDPhysicsAllParameters::fromJson(nlohmann::json& physicsJsonParams)
{
	physicsParams.fromJson(physicsJsonParams["PhysicsParams"]);
	collisionParams.fromJson(physicsJsonParams["CollisionParams"]);
	if (collisionParams.allowVolumetricCollision) {
		collisionParams.volCollisionParams.fromJson(physicsJsonParams["VolCollisionParams"]);
	}
	return true;
}

bool GAIA::PBDPhysicsAllParameters::toJson(nlohmann::json& physicsJsonParams)
{
	physicsParams.toJson(physicsJsonParams["PhysicsParams"]);
	collisionParams.toJson(physicsJsonParams["CollisionParams"]); 
	if (collisionParams.allowVolumetricCollision) {
		collisionParams.volCollisionParams.toJson(physicsJsonParams["VolCollisionParams"]);
	}
	return true;
}

#include <CollisionDetector/CollisionGeometry.h>
#include <CollisionDetector/TriMeshCollisionGeometry.h>

#include "VBDClothPhysics.h"
#include <Timer/Timer.h>
#include <IO/FileIO.h>
#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>
// #define TURN_ON_SINGLE_THREAD_DEBUG
#include <Parallelization/CPUParallelization.h>

#include "VBDClothDeformer.h"
#include <limits>

#include <CollisionDetector/TriangleTriangleIntersection.h>

using namespace GAIA;

void GAIA::VBDClothSimulationFramework::initialize()
{
	BaseClothPhsicsFramework::initialize();

	orgPosistions.resize(numClothes());

	if (physicsParams().useNewton)
	{
		initializeNewton();
	}
	else {
		initializeParallelGroups();
	}

	// load deformers
	auto deformerParams = physicsJsonParams["Deformers"];

	for (auto deformerParam : deformerParams) {
		deformers.push_back(loadClothDeformers(*this, deformerParam));
	}

	if (pViewerParams->enableViewer)
	{
		initializeViewer();
	}

	loadConstraints();

	// initialize collision data
	clothClothVFContactQueyResults.resize(triMeshesAll.size());
	clothClothEEContactQueyResults.resize(baseTriMeshesForSimulation.size());
	clothClothFVContactQueyResults.resize(baseTriMeshesForSimulation.size());

	for (IdType iMesh = 0; iMesh < triMeshesAll.size(); iMesh++) {
		clothClothVFContactQueyResults[iMesh].resize(triMeshesAll[iMesh]->numVertices());
	}
	for (IdType iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		clothClothEEContactQueyResults[iMesh].resize(baseTriMeshesForSimulation[iMesh]->numEdges());

		// not used for now
		//clothClothFVContactQueyResults.resize(baseTriMeshesForSimulation[iMesh]->numFaces());
	}


}

void GAIA::VBDClothSimulationFramework::initializeNewton()
{
	pNewtonAssembler = std::make_shared<TriMeshNewtonAssembler>();
	pNewtonAssembler->initialize(baseTriMeshesForSimulation, physicsParams().newtonSolverType);
}

void GAIA::VBDClothSimulationFramework::initializeParallelGroups()
{
	// generate parallelization groups
	size_t numberOfParallelGroups = 0;
	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		if (numberOfParallelGroups < baseTriMeshesForSimulation[iMesh]->verticesColoringCategories().size())
		{
			numberOfParallelGroups = baseTriMeshesForSimulation[iMesh]->verticesColoringCategories().size();
		}
	}

	vertexParallelGroups.resize(numberOfParallelGroups);

	for (int iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		auto pMesh = baseTriMeshesForSimulation[iMesh];
		size_t numColors = pMesh->verticesColoringCategories().size();
		std::vector<IdType> availableColors(numColors, -1);

		for (int iColor = 0; iColor < numColors; iColor++)
		{
			availableColors[iColor] = iColor;
		}

		for (int iColor = 0; iColor < numColors; iColor++)
		{
			const std::vector<IdType>& currentColorGroup = pMesh->verticesColoringCategories()[iColor];

			int smallestGroupId = findSmallestParallelGroup(availableColors, vertexParallelGroups, false);
			// add the current color group to this parallel group
			for (int iVertex = 0; iVertex < currentColorGroup.size(); iVertex++)
			{
				vertexParallelGroups[smallestGroupId].push_back(iMesh);
				vertexParallelGroups[smallestGroupId].push_back(currentColorGroup[iVertex]);
				pMesh->globalColors[currentColorGroup[iVertex]] = smallestGroupId;
			}
			// color group from the same mesh cannot be added to this group again
			// std::cout << "Color " << iColor << " of mesh " << iMesh << " was added to parallel group " << smallestGroupId << "\n";
		}
	}
}

IdType GAIA::VBDClothSimulationFramework::findSmallestParallelGroup(std::vector<IdType>& availableColors, const std::vector<std::vector<IdType>>& vertexParallelGroups,
	bool removeSmallestAvailableColor) {
	assert(vertexParallelGroups.size());
	assert(availableColors.size());

	size_t smallestParallelGroupSize = vertexParallelGroups[availableColors[0]].size();
	int smallestParallelGroupId = availableColors[0];
	IdType smallestI = 0;
	for (IdType i = 1; i < availableColors.size(); i++)
	{
		if (vertexParallelGroups[availableColors[i]].size() < smallestParallelGroupSize)
		{
			smallestParallelGroupSize = vertexParallelGroups[availableColors[i]].size();
			smallestParallelGroupId = availableColors[i];
			smallestI = i;
		}

	}
	if (removeSmallestAvailableColor)
	{
		availableColors.erase(availableColors.begin() + smallestI);
	}
	return smallestParallelGroupId;
}

void GAIA::VBDClothSimulationFramework::parseRunningParameters()
{
	BaseClothPhsicsFramework::parseRunningParameters();

	pTargetMatcherParameter = std::make_shared<MeshClosestPointQueryParameters>();
	pTargetMatcherParameter->fromJson(physicsJsonParams["MeshClosestPointQueryParams"]);
}

void GAIA::VBDClothSimulationFramework::loadConstraints()
{
	// load vertex constraints

	if (physicsJsonParams["PhysicsParams"].contains("VertexFaceAttachmentConstraintParams"))
	{
		VertexFaceAttachmentConstraintParams vfAttachmentConstraintParams;
		vfAttachmentConstraintParams.fromJson(physicsJsonParams["PhysicsParams"]["VertexFaceAttachmentConstraintParams"]);

		RELEASE_ASSERT(vfAttachmentConstraintParams.vertexConstraintVertexInfos.size() % 2 == 0);
		RELEASE_ASSERT(vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos.size() % 4 == 0);
		size_t numConstraints = vfAttachmentConstraintParams.vertexConstraintVertexInfos.size() / 2;

		for (size_t iConstraint = 0; iConstraint < numConstraints; iConstraint++)
		{
			IdType meshIdVSide = vfAttachmentConstraintParams.vertexConstraintVertexInfos[iConstraint * 2];
			IdType vertexId = vfAttachmentConstraintParams.vertexConstraintVertexInfos[iConstraint * 2 + 1];

			RELEASE_ASSERT(meshIdVSide < numSimulationMeshes());
			VBDBaseTriMesh* pMeshVSide = &getSimulatedMesh(meshIdVSide);
			pMeshVSide->vertexConstraints[vertexId].attachmentConstraints.emplace_back();

			VertexFaceAttachmentConstraint& vertexAttachmentConstraint = pMeshVSide->vertexConstraints[vertexId].attachmentConstraints.back();

			IdType targetMeshId = vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos[iConstraint * 4];

			// targetMeshId = -1 means attaching to a fixed point
			vertexAttachmentConstraint.vertexOrder = 3;
			vertexAttachmentConstraint.attachedMeshFaceSide = targetMeshId;
			vertexAttachmentConstraint.attachedMeshVSdie = meshIdVSide;
			vertexAttachmentConstraint.vertexId = vertexId;

			if (vfAttachmentConstraintParams.vertexConstraintStiffness.size())
			{
				vertexAttachmentConstraint.constraintStiffness = vfAttachmentConstraintParams.vertexConstraintStiffness[iConstraint];
			}

			Vec3I face = {
				vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos[iConstraint * 4 + 1],
				vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos[iConstraint * 4 + 2],
				vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos[iConstraint * 4 + 3]
			};
			vertexAttachmentConstraint.attachedFaceVIds = face;

			Vec3 attachedPos = {
				vfAttachmentConstraintParams.vertexConstraintBarys[iConstraint][0],
				vfAttachmentConstraintParams.vertexConstraintBarys[iConstraint][1],
				vfAttachmentConstraintParams.vertexConstraintBarys[iConstraint][2]
			};
			vertexAttachmentConstraint.attachedPos = attachedPos;

			if (targetMeshId < numSimulationMeshes() && targetMeshId >= 0)
				// another side is a simulated mesh, we need to add constriants to the other side too
			{
				VBDBaseTriMesh* pMeshFSide = &getSimulatedMesh(targetMeshId);
				for (IdType iFV = 0; iFV < 3; iFV++)
				{
					if (attachedPos[iFV] == 0.f)
					{
						continue;
					}
					IdType faceVId = vfAttachmentConstraintParams.vertexConstraintTargetFacesInfos[iConstraint * 4 + 1 + iFV];
					pMeshFSide->vertexConstraints[faceVId].attachmentConstraints.emplace_back();
					VertexFaceAttachmentConstraint& vertexAttachmentConstraintFSide =
						pMeshFSide->vertexConstraints[faceVId].attachmentConstraints.back();
					vertexAttachmentConstraintFSide = vertexAttachmentConstraint;
					vertexAttachmentConstraintFSide.vertexOrder = iFV;
					vertexAttachmentConstraintFSide.attachedFaceVIds = face;
					vertexAttachmentConstraintFSide.attachedPos = attachedPos;

				}
			}
		}
	}


}

TriMeshFEM::SharedPtr VBDClothSimulationFramework::initializeMaterial(ObjectParams::SharedPtr objParam, BasePhysicsParams::SharedPtr physicsParaemters, BaseClothPhsicsFramework* pPhysics)
{
	TriMeshFEM::SharedPtr pMesh;
	if (objParam->materialName == "StVK_triMesh")
	{
		VBDTriMeshStVk::SharedPtr pMeshStVK = std::make_shared<VBDTriMeshStVk>();
		pMeshStVK->initialize(std::static_pointer_cast<TriMeshParams>(objParam), this);
		pMesh = pMeshStVK;
	}

	return pMesh;
}

std::shared_ptr<ObjectParamsList> GAIA::VBDClothSimulationFramework::createObjectParamsList()
{
	return std::make_shared<ObjectParamsListVBDCloth>();
}

std::shared_ptr<BasePhysicsParams> GAIA::VBDClothSimulationFramework::createPhysicsParams()
{
	return std::make_shared<VBDClothPhysicsParameters>();
}

std::shared_ptr<RunningTimeStatistics> GAIA::VBDClothSimulationFramework::createRunningTimeStatistics()
{
	return std::make_shared<VBDClothRuntimeStatistics>();
}

void GAIA::VBDClothSimulationFramework::runStepSequential()
{
	applyInitialGuess();
	ConvergenceStats convergenceStat(baseTriMeshesForSimulation.size());

	iIter = -1;
	FloatingType energy = evaluateEnergy(convergenceStat);
	printStats(convergenceStat, energy);
	FloatingType energyPrev = energy;
	FloatingType avgGradPrev = convergenceStat.avgGradNorm_allClothes;
	FloatingType stepSize = physicsParams().stepSize;

	int stepsToLastLineSearch = 0;

	convergenceStats.clear();
	convergenceStats.push_back(convergenceStat);
	for (iIter = 0; iIter < physicsParams().maxNumIteration; iIter++)
	{
		for (size_t iCloth = 0; iCloth < numClothes(); iCloth++)
		{
			VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
			// mesh.solverIteration();
		}

		//if (!(iIter % physicsParams().evaluationSteps))
		//{

		//	clearGradient();
		//	accumulateMaterialGradient();

		//	energy = evaluateEnergy(convergenceStat);

		//	convergenceStats.push_back(convergenceStat);
		//	printStats(convergenceStat, energy);

		//	for (size_t iCloth = 0; iCloth < numClothes(); iCloth++)
		//	{
		//		VBDTriMeshStVk& mesh = getSimulatedMesh(iCloth);

		//	}

		//	if (determineConvergence(convergenceStat, energyPrev, avgGradPrev, stepSize))
		//	{
		//		std::cout << "Converged at iter: " << iIter << std::endl;
		//		break;
		//	}
		//	avgGradPrev = convergenceStat.avgGradNorm_allClothes;
		//	energyPrev = energy;
		//}
	}

	if (physicsParams().outputStatistics)
	{
		std::string statisticsOutOutPath = outputFolder + "/Statistics/" + "/Frame_" + MF::STR::padToNew(std::to_string(frameId), 4, '0') + ".json";
		saveStats(statisticsOutOutPath, convergenceStats);
	}

	updateVelocity();
}

void GAIA::VBDClothSimulationFramework::runStep()
{
	//runStep_CDEveryColor();
	if (physicsParams().useNewton)
	{
		runStep_Newton();
	}
	else {
		runStep_CDEveryIter();
	}
}

void GAIA::VBDClothSimulationFramework::runStep_CDEveryColor()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});


		applyInitialGuess();
		applyDeformers();
		TICK(timeCsmpMaterialSolve);
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
				std::cout << "iIter: " << iIter << std::endl;
				});
			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
			{
				//pClothContactDetector->updateBVH();
				updateBVH(iGroup);

				const std::vector<IdType>& parallelGroup = vertexParallelGroups[iGroup];

				size_t numVertices = parallelGroup.size() / 2;
				cpu_parallel_for(0, numVertices, [&](int iV) {
					//for (size_t iV = 0; iV < numVertices; iV++) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];
					auto pMesh = &getSimulatedMesh(iMesh);

					if (!pMesh->fixedMask[vId])
					{
						//pMesh->VBDStep(vId);
						VBDStepWithCollisionDetection(iMesh, vId, true);
					}
					});


				// copy new positions to the current positions to ensure GS iteration
				cpu_parallel_for(0, numVertices, [&](int iV) {
					//for (size_t iV = 0; iV < numVertices; iV++) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];
					auto pMesh = &getSimulatedMesh(iMesh);
					if (!pMesh->fixedMask[vId])
					{
						pMesh->vertex(vId) = pMesh->vertexPositionNext(vId);
						//std::cout << "pMesh->vertex(" << vId << "):" << pMesh->vertexPositionNext(vId).transpose() << std::endl;
					}
					});
			}
		} // iteration
		TOCK_STRUCT(runtimeStatistics(), timeCsmpMaterialSolve);

		// detectIntersections();
		//applyAcceleration(1.f);

		updateVelocity();
	} // substep
}

void GAIA::VBDClothSimulationFramework::runStep_CDEveryIter()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});

		ChebysevAccelerator accelerator(physicsParams().acceleratorRho);
		FloatingType omega = 1.f;
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{

			debugOperation(DEBUG_LVL_DEBUG, [&]() {
				std::cout << "iIter: " << iIter << std::endl;
				});
			TICK(timeCsmpUpdatingCollisionInfoDCD);

			if (iIter == 0)
			{
				applyDeformers();
				initializeActiveCollisionMask();
			}
			updateCollider();

			// Do collision detection if colliders are updated
			for (const auto& colliderMesh : colliderTriMeshes) {
				collisionDetectionRequired = (colliderMesh->updated || collisionDetectionRequired);
			}
			if (collisionDetectionRequired || iIter == 0)
				// collision detection is required when one vertex moves out of their conservative bounds
			{
				if (physicsParams().handleCollision)
				{
					applyCollisionDetection();
					runtimeStatistics().numberCollisionDetection++;
				}
				else
				{
					for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
					{
						VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
						pMesh->vertexConvervativeBounds.setConstant(std::numeric_limits<FloatingType>::max());
					}
				}
			}
			TOCK_STRUCT(runtimeStatistics(), timeCsmpUpdatingCollisionInfoDCD);

			if (physicsParams().applyAcceleration)
			{
				for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
				{
					VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
					pMesh->positionPrevIter = pMesh->positions();
				}
			}

			TICK(timeCsmpInitialStep);
			if (iIter == 0)
			{
				applyInitialGuess();
				// record initial residuals
				if (physicsParams().evaluateConvergence)
				{
					evaluateConvergence();
				}
			}
			TOCK_STRUCT(runtimeStatistics(), timeCsmpInitialStep);

			TICK(timeCsmpMaterialSolve);
			setSimulatedMeshToUpToDateStatus(false);

			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
			{
				//pClothContactDetector->updateBVH();

				const std::vector<IdType>& parallelGroup = vertexParallelGroups[iGroup];

				size_t numVertices = parallelGroup.size() / 2;
				cpu_parallel_for(0, numVertices, [&](int iV) {
					//for (size_t iV = 0; iV < numVertices; iV++) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];
					auto pMesh = &getSimulatedMesh(iMesh);

					if (!pMesh->fixedMask[vId])
					{
						//pMesh->VBDStep(vId);
						VBDStepWithExistingCollisions(iMesh, vId);
					}
					});

				// copy new positions to the current positions to ensure GS iteration
				cpu_parallel_for(0, numVertices, [&](int iV) {
					//for (size_t iV = 0; iV < numVertices; iV++) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];
					auto pMesh = &getSimulatedMesh(iMesh);
					if (!pMesh->fixedMask[vId])
					{
						pMesh->vertex(vId) = pMesh->vertexPositionNext(vId);
						//std::cout << "pMesh->vertex(" << vId << "):" << pMesh->vertexPositionNext(vId).transpose() << std::endl;
					}
					});
			}
			setSimulatedMeshToUpToDateStatus(true);

			if (physicsParams().applyAcceleration)
			{
				omega = accelerator.getAcceleratorOmega(iIter+1, omega);
				applyAcceleration(omega);
			}
			TOCK_STRUCT(runtimeStatistics(), timeCsmpMaterialSolve);

			if (physicsParams().evaluateConvergence)
			{
				evaluateConvergence();
			}
			if (physicsParams().saveIntermediateResults)
			{
				saveIntermediateResults();
			}
		} // iteration

		// detectIntersections();
		
		updateVelocity();
	} // substep
}

void GAIA::VBDClothSimulationFramework::runStep_Newton()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});


		bool newtonHessianPatternChanged = false;
		NFloatingType newtonEnergy{};

		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{

			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
				std::cout << "iIter: " << iIter << std::endl;
				});

			TICK(timeCsmpUpdatingCollisionInfoDCD);
			if (physicsParams().handleCollision) {
				if (collisionDetectionRequired || iIter == 0)
					// collision detection is required when one vertex moves out of their conservative bounds
				{
					std::cout << "Collision detection required at iter: " << iIter << std::endl;
					applyCollisionDetection();
					pNewtonAssembler->analyzeCollision(clothClothVFContactQueyResults, clothClothEEContactQueyResults);
					newtonHessianPatternChanged = true;
					runtimeStatistics().numberCollisionDetection++;
				}
			}
			TOCK_STRUCT(runtimeStatistics(), timeCsmpUpdatingCollisionInfoDCD);

			TICK(timeCsmpInitialStep);
			if (iIter == 0)
			{
				applyInitialGuess();
				applyDeformers();
				updateCollider();

				if (physicsParams().evaluateConvergence)
				{
					evaluateConvergence();
				}
			}
			TOCK_STRUCT(runtimeStatistics(), timeCsmpInitialStep);

			setSimulatedMeshToUpToDateStatus(false);

			TICK(timeCsmpMaterialSolve);
			fillNewtonSystem();
			if (newtonHessianPatternChanged)
			{
				newtonHessianPatternChanged = false;
				pNewtonAssembler->solve(true, physicsParams().handleCollision);
			}
			else
			{
				pNewtonAssembler->solve(false, physicsParams().handleCollision);
			}
			newtonDx = pNewtonAssembler->Ndx.cast<FloatingType>();

			if (physicsParams().useLineSearch) {
				if (iIter == 0) {
					NFloatingType eInertia{};
					NFloatingType eElastic{};
					NFloatingType eBending{};
					NFloatingType eContact{};
					// Elastic is already computed during force and hessian evaluation, set elasticReady to true
					newtonEnergy = newtonEvaluateMeritEnergy(eInertia, eElastic, eBending, eContact, true);
					debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
						std::cout << "energy: " << newtonEnergy << ", inertia: " << eInertia << ", elastic: " << eElastic << ", bending: " << eBending << ", contact: " << eContact << std::endl;
						});
				}
				FloatingType stepSizeNew = 1.f;
				newtonEnergy = newtonLineSearch(newtonDx, newtonEnergy, stepSizeNew, physicsParams().c,
					physicsParams().tau, physicsParams().backtracingLineSearchMaxIters, stepSizeNew);

			}
			else {
				pNewtonAssembler->updatePositions(newtonDx);
			}
			// TODO: move culling before line search
			if (physicsParams().handleCollision) {
				newtonConservativeStepCulling();
			}

			setSimulatedMeshToUpToDateStatus(true);
			TOCK_STRUCT(runtimeStatistics(), timeCsmpMaterialSolve);
			if (physicsParams().evaluateConvergence)
			{
				evaluateConvergence();
			}
		} // iteration

		// detectIntersections();
		//applyAcceleration(1.f);

		updateVelocity();
	} // substep
}

//void GAIA::VBDClothSimulationFramework::runStep_blockJacobi()
//{
//	applyInitialGuess();
//	ConvergenceStats convergenceStat(baseTriMeshesForSimulation.size());
//
//	iIter = -1;
//	FloatingType energy = evaluateEnergy(convergenceStat);
//	printStats(convergenceStat, energy);
//	FloatingType energyPrev = energy;
//	FloatingType avgGradPrev = convergenceStat.avgGradNorm_allClothes;
//	FloatingType stepSize = physicsParams().stepSize;
//
//	int stepsToLastLineSearch = 0;
//
//	convergenceStats.clear();
//	for (iIter = 0; iIter < physicsParams().maxNumIteration; iIter++)
//	{
//		clearGradient();
//		accumulateMaterialGradient();
//		accumulateRegistrationGradient();
//
//		if (!(iIter % physicsParams().contactDetectionIters) && physicsParams().handleCollision)
//		{
//			collisionDetection();
//		}
//		accumulateCollisionGradient();
//
//		preconditionGradients();
//
//		if (iIter < physicsParams().numLineSearchSteps && stepSize > physicsParams().minStepSizeGD)
//		{
//			stepSize = backTracingLineSearchVBD(energy, stepSize, physicsParams().c,
//				physicsParams().tau, physicsParams().backtracingLineSearchMaxIters, physicsParams().minStepSizeGD, convergenceStat);
//
//			std::cout << " - Applied initial line search at iter: " << iIter << ", step size changed to: " << stepSize << std::endl;
//		}
//		else {
//			if (stepsToLastLineSearch > physicsParams().lineSearchEveryNumSteps)
//			{
//				stepSize = backTracingLineSearchVBD(energy, physicsParams().stepSize, physicsParams().c,
//					physicsParams().tau, physicsParams().backtracingLineSearchMaxIters, physicsParams().minStepSizeGD, convergenceStat);
//				stepsToLastLineSearch = 0;
//			}
//			else
//			{
//				applyGradientDescent(stepSize);
//				++stepsToLastLineSearch;
//			}
//		}
//
//		if (physicsParams().saveIntermediateResults)
//		{
//			saveIntermediateResults();
//		}
//
//		if (!(iIter % physicsParams().evaluationSteps)
//			|| (iIter < physicsParams().numLineSearchSteps))
//		{
//			energy = evaluateEnergy(convergenceStat);
//
//			//if ((energyPrev < energy 
//			//	|| ) 
//			//	&& stepSize > physicsParams().minStepSizeGD)
//			//{
//			//	// revert previous step
//			//	revertGradientDescent(stepSize);
//
//			//	stepSize = backTracingLineSearchVBD(energy, stepSize, physicsParams().c,
//			//		physicsParams().tau, physicsParams().backtracingLineSearchMaxIters, physicsParams().minStepSizeGD);
//
//			//	std::cout << " - Applied line search at iter: " << iIter << ", step size changed to: " << stepSize << std::endl;
//			//	stepsToLastLineSearch = 0;
//			//}
//
//			convergenceStats.push_back(convergenceStat);
//			printStats(convergenceStat, energy);
//
//			if (determineConvergence(convergenceStat, energyPrev, avgGradPrev, stepSize))
//			{
//				std::cout << "Converged at iter: " << iIter << std::endl;
//				break;
//			}
//			avgGradPrev = convergenceStat.avgGradNorm_allClothes;
//			energyPrev = energy;
//		}
//	}
//
//	if (physicsParams().outputStatistics)
//	{
//		std::string statisticsOutOutPath = outputFolder + "/Statistics/" + "/Frame_" + MF::STR::padToNew(std::to_string(frameId), 4, '0') + ".json";
//		saveStats(statisticsOutOutPath, convergenceStats);
//	}
//
//	updateVelocity();
//}

void GAIA::VBDClothSimulationFramework::updateBVH(int groupId)
{
	if (frameId
		&& (frameId % physicsParams().contactBVHReconstructionSteps == 0)
		&& (iIter == 0)
		&& (groupId == 0)
		)
	{
		pClothContactDetector->updateBVH(RTC_BUILD_QUALITY_LOW);
	}
	else
	{
		pClothContactDetector->updateBVH(RTC_BUILD_QUALITY_REFIT);
	}
}

void GAIA::VBDClothSimulationFramework::applyCollisionDetection()
{
	TICK(timeCsmpUpdatingBVHDCD);
	updateBVH();
	pClothContactDetector->resetFaceContactInfo();
	TOCK_STRUCT(runtimeStatistics(), timeCsmpUpdatingBVHDCD);

	TICK(timeCsmpColDetectDCD);
	for (size_t iMesh = 0; iMesh < triMeshesAll.size(); iMesh++)
	{
		auto* pMesh = &getMesh(iMesh);
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {
			pClothContactDetector->contactQueryVF(iMesh, iV, &getVFContactResult(iMesh, iV));
			});
	}
	// only do ee collision detection if the mesh is simulated
	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		auto* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numEdges(), [&](int iE) {
			pClothContactDetector->contactQueryEE(iMesh, iE, &getEEContactResult(iMesh, iE));
			});

		//cpu_parallel_for(0, pMesh->numFaces(), [&](int fId) {

		//	ClothVFContactQueryResult& contactResultFV_ = getFVContactResult(iMesh, fId);
		//	pClothContactDetector->contactQueryFV(iMesh, fId, &contactResultFV_, -1);
		//	// The mesh has been modified, therefore they aren't suppose to match
		//	// But at least they should match in the first parallelization group
		//	CFloatingType faceMinDisToVertices = pClothContactDetector->getFaceMinDis(iMesh, fId);
		//	if ((contactResultFV_.minDisToPrimitives != faceMinDisToVertices)
		//		&& contactResultFV_.minDisToPrimitives != std::numeric_limits<FloatingType>::max()
		//		) {
		//		std::cout << "fId: " << fId << " | contactResultFV_.minDisToPrimitives: " << contactResultFV_.minDisToPrimitives
		//			<< " | faceMinDisToVertices: " << faceMinDisToVertices << std::endl;

		//		assert(contactResultFV_.minDisToPrimitives >= faceMinDisToVertices);
		//	}
		//});
	}
	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		auto* pMesh = &getSimulatedMesh(iMesh);
		// recomputed the conservative bounds
		cpu_parallel_for(0, pMesh->numVertices(), [&](int vId) {
			FloatingType dMin = pClothContactDetectorParameters->maxQueryDis;
			ClothVFContactQueryResult& contactResult = getVFContactResult(iMesh, vId);
			dMin = std::min(dMin, contactResult.minDisToPrimitives);
			
			bool activelyColliding = contactResult.contactPts.size() > 0;

			for (size_t iNeiF = 0; iNeiF < pMesh->numNeiFaces(vId); iNeiF++)
			{
				IdType fId = pMesh->getVertexIthNeiFace(vId, iNeiF);
				CFloatingType faceMinDisToVertices = pClothContactDetector->getFaceMinDis(iMesh, fId);
				dMin = std::min(dMin, faceMinDisToVertices);
				activelyColliding = activelyColliding || pClothContactDetector->getFaceContactInfo(iMesh, fId).size() > 0;
			}

			for (size_t iNeiE = 0; iNeiE < pMesh->numNeiEdges(vId); iNeiE++)
			{
				IdType eId = pMesh->getVertexIthNeiEdge(vId, iNeiE);
				ClothEEContactQueryResult& contactResultEE = getEEContactResult(iMesh, eId);
				dMin = std::min(dMin, contactResultEE.minDisToPrimitives);
				activelyColliding = activelyColliding || contactResultEE.contactPts.size() > 0;
			}

			//RELEASE_ASSERT_WITH_FUNC(dMin - physicsParams().thickness * 0.5f >= 0, 
			//	[&]() {
			//		std::cout << "Vertex Id: " << vId << '\n'
			//			<< "dMin: " << dMin << " | physicsParams().thickness * 0.5f: " << physicsParams().thickness * 0.5f << std::endl;
			//	});
			dMin = std::max(0.f, dMin - physicsParams().thickness);

			pMesh->vertexConvervativeBounds[vId] = dMin * physicsParams().conservativeStepRelaxation;
			
			// activeCollisionMask can only be set to true, but not false
			if (activelyColliding)
			{
				pMesh->activeCollisionMask(vId) = true;
			}
			});

		pMesh->positionsAtPrevCollisionDetection = pMesh->positions();
	}
	TOCK_STRUCT(runtimeStatistics(), timeCsmpColDetectDCD);

	collisionDetectionRequired = false;
	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
		std::vector<FloatingType> vfdistances{};
		std::vector<FloatingType> self_vfdistances{};
		std::vector<std::array<int, 7>> vfpairs{};
		for (int iV = 0; iV < clothClothVFContactQueyResults[0].size(); iV++)
		{
			for (int iC = 0; iC < clothClothVFContactQueyResults[0][iV].numContactPoints(); iC++)
			{
				vfdistances.push_back(clothClothVFContactQueyResults[0][iV].contactPts[iC].d);
				const auto vid = clothClothVFContactQueyResults[0][iV].contactPts[iC].contactVertexId;
				const auto fid = clothClothVFContactQueyResults[0][iV].contactPts[iC].contactFaceId;
				const auto mesh_id1 = clothClothVFContactQueyResults[0][iV].contactPts[iC].contactVertexSideMeshId;
				const auto mesh_id2 = clothClothVFContactQueyResults[0][iV].contactPts[iC].contactFaceSideMeshId;
				const auto& mesh2 = getMesh(mesh_id2);
				const int fvi[3] = { mesh2.facePosVId(fid,0), mesh2.facePosVId(fid,1), mesh2.facePosVId(fid,2) };
				vfpairs.push_back({ mesh_id1, vid, mesh_id2, fid, fvi[0], fvi[1], fvi[2] });
				if (mesh_id1 == 0 && mesh_id2 == 0) {
					self_vfdistances.push_back(clothClothVFContactQueyResults[0][iV].contactPts[iC].d);
				}
			}
		}
		std::vector<FloatingType> eedistances{};
		std::vector<FloatingType> self_eedistances{};
		std::vector<std::array<int, 8>> eepairs{};
		for (int iE = 0; iE < clothClothEEContactQueyResults[0].size(); iE++)
		{
			for (int iC = 0; iC < clothClothEEContactQueyResults[0][iE].numContactPoints(); iC++)
			{
				eedistances.push_back(clothClothEEContactQueyResults[0][iE].contactPts[iC].d);
				const auto eid1 = clothClothEEContactQueyResults[0][iE].contactPts[iC].contactEdgeId1;
				const auto eid2 = clothClothEEContactQueyResults[0][iE].contactPts[iC].contactEdgeId2;
				const auto mesh_id1 = clothClothEEContactQueyResults[0][iE].contactPts[iC].contactMeshId1;
				const auto mesh_id2 = clothClothEEContactQueyResults[0][iE].contactPts[iC].contactMeshId2;
				const auto& mesh1 = getMesh(mesh_id1);
				const auto& mesh2 = getMesh(mesh_id2);
				const int evi1[2] = { mesh1.getEdgeInfo(eid1).eV1, mesh1.getEdgeInfo(eid1).eV2 };
				const int evi2[2] = { mesh2.getEdgeInfo(eid2).eV1, mesh2.getEdgeInfo(eid2).eV2 };
				eepairs.push_back({ mesh_id1, eid1, evi1[0], evi1[1], mesh_id2, eid2, evi2[0], evi2[1] });
				if (mesh_id1 == 0 && mesh_id2 == 0) {
					self_eedistances.push_back(clothClothEEContactQueyResults[0][iE].contactPts[iC].d);
				}
			}
		}
		const FloatingType min_vfdistance = vfdistances.size() > 0 ? *std::min_element(vfdistances.begin(), vfdistances.end()) : -1;
		const FloatingType min_eedistance = eedistances.size() > 0 ? *std::min_element(eedistances.begin(), eedistances.end()) : -1;
		const FloatingType min_self_vfdistance = self_vfdistances.size() > 0 ? *std::min_element(self_vfdistances.begin(), self_vfdistances.end()) : -1;
		const FloatingType min_self_eedistance = self_eedistances.size() > 0 ? *std::min_element(self_eedistances.begin(), self_eedistances.end()) : -1;
		std::cout << "min_vfdistance: " << min_vfdistance << " | min_eedistance: " << min_eedistance << std::endl;
		std::cout << "min_self_vfdistance: " << min_self_vfdistance << " | min_self_eedistance: " << min_self_eedistance << std::endl;
		nlohmann::json j;
		j["vfdistances"] = vfdistances;
		j["eedistances"] = eedistances;
		j["vfpairs"] = vfpairs;
		j["eepairs"] = eepairs;
		std::string outPath = getDebugFolder() + "/dist_" + std::to_string(frameId) + "_" + std::to_string(substep) + "_" + std::to_string(iIter) + ".json";
		std::ofstream o(outPath);
		o << j << std::endl;
		o.close();
		});
}

void GAIA::VBDClothSimulationFramework::initializeActiveCollisionMask()
{
	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		pMesh->activeCollisionMask.setConstant(false);
	}
}

void GAIA::VBDClothSimulationFramework::applyAcceleration(CFloatingType omega)
{
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iMesh);
		// mesh.positionsPrev = mesh.positions();
		if (omega > 1.f)
		{
			cpu_parallel_for(0, mesh.numVertices(), [&](int iV) {

				if (!mesh.activeCollisionMask(iV))
				{
					Vec3 posNew = omega * (mesh.vertex(iV) - mesh.positionPrevIter.col(iV)) + mesh.positionPrevIter.col(iV);
					bool truncationNeeded = applyConservativeBoundTruncation(iMesh, iV, posNew);
					collisionDetectionRequired = truncationNeeded;
					mesh.vertex(iV) = posNew;
				}

			});
		}
		mesh.positionPrevPrevIter = mesh.positionPrevIter;
	}
}

void GAIA::VBDClothSimulationFramework::VBDStepWithCollisionDetection(IdType iMesh, IdType vId, bool updateCollisions)
{
	VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);

	Mat3 hessian = Mat3::Zero();
	Vec3 force = Vec3::Zero();
	if (pMesh->fixedMask[vId])
	{
		return;
	}

	FloatingType dMin = pClothContactDetectorParameters->maxQueryDis;

	// vertex collision detection
	ClothVFContactQueryResult& contactResult = getVFContactResult(iMesh, vId);
	if (updateCollisions)
	{
		pClothContactDetector->contactQueryVF(iMesh, vId, &contactResult);
	}

	for (size_t iContact = 0; iContact < contactResult.numContactPoints(); iContact++)
	{
		// if collision detection is not applied update the contact point info
		// vertex is from the v side, so the vertex order is always 3
		CFloatingType dis = accumulateVFContactForceAndHessian(iContact, 3, &contactResult, !updateCollisions, true, force, hessian);
	}
	dMin = std::min(dMin, contactResult.minDisToPrimitives);

	// neighbor faces collision detection
	for (size_t iNeiF = 0; iNeiF < pMesh->numNeiFaces(vId); iNeiF++)
	{
		IdType fId = pMesh->getVertexIthNeiFace(vId, iNeiF);
		ClothVFContactQueryResult& contactResultFV = getFVContactResult(iMesh, fId);

		if (updateCollisions)
		{
			pClothContactDetector->contactQueryFV(iMesh, fId, &getFVContactResult(iMesh, fId), vId);
		}

		dMin = std::min(dMin, contactResultFV.minDisToPrimitives);

		for (size_t iContact = 0; iContact < contactResultFV.numContactPoints(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			CFloatingType dis = accumulateVFContactForceAndHessian(iContact, pMesh->getVertexIthNeiFaceOrder(vId, iNeiF), &contactResultFV,
				!updateCollisions, true, force, hessian);
			//dMin = std::min(dis, dMin);
		}
	}

	// neighbor edges collision detection

	for (size_t iNeiE = 0; iNeiE < pMesh->numNeiEdges(vId); iNeiE++)
	{
		IdType eId = pMesh->getVertexIthNeiEdge(vId, iNeiE);
		ClothEEContactQueryResult& contactResultEE = getEEContactResult(iMesh, eId);
		if (updateCollisions)
		{
			CFloatingType dis = pClothContactDetector->contactQueryEE(iMesh, eId, &contactResultEE);
		}
		for (size_t iContact = 0; iContact < contactResultEE.numContactPoints(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			CFloatingType dis = accumulateEEContactForceAndHessian(iContact, pMesh->getVertexIthNeiEdgeOrder(vId, iNeiE), &contactResultEE,
				!updateCollisions, true, force, hessian);
			// dMin = std::min(dis, dMin);
		}

		dMin = std::min(dMin, contactResultEE.minDisToPrimitives);

	}

	// accumulate the inertia and elasticity force and hessian
	pMesh->accumlateInertiaForceAndHessian(vId, force, hessian);
	pMesh->accumlateMaterialForceAndHessian(vId, force, hessian);
	// accumulate the constraints' force and hessian
	accumulateConstraintForceAndHessian(iMesh, vId, force, hessian);

	// add external force
	force += pMesh->vertexExternalForces.col(vId);

	// Solve linear system
	Vec3 descentDirection;
	FloatingType stepSize = physicsParams().stepSize;
	bool solverSuccess = CuMatrix::solve3x3_psd_stable(hessian.data(), force.data(), descentDirection.data());

	if (!solverSuccess)
	{
		std::cout << "Solver failed at vertex " << vId << std::endl;
	}
	else
	{
		Vec3 descentStep = stepSize * descentDirection;
		FloatingType descentStepSize = descentStep.norm();
		if (descentStepSize > physicsParams().conservativeStepRelaxation * dMin) {
			FloatingType descentStepSizeNew = physicsParams().conservativeStepRelaxation * dMin;
			descentStep = descentStep * descentStepSizeNew / descentStepSize;
			// std::cout << "Step truncation triggered for vertex " << vId << ", from " << descentStepSize 
			// 	<< " to " << descentStepSizeNew 
			// 	<< " at frame" << frameId << " iter " << iIter 
			// 	<< "\n";
		}

		// update vertex position to PositionNext to ensure no racing condition
		pMesh->vertexPositionNext(vId) = pMesh->vertex(vId) + descentStep;
		//pMesh->vertex(vId) = pMesh->vertex(vId) + stepSize * descentDirection;
	}
}

void GAIA::VBDClothSimulationFramework::VBDStepWithExistingCollisions(IdType iMesh, IdType vId)
{
	VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);

	Mat3 hessian = Mat3::Zero();
	Vec3 force = Vec3::Zero();

	ClothVFContactQueryResult& contactResult = getVFContactResult(iMesh, vId);

	for (size_t iContact = 0; iContact < contactResult.numContactPoints(); iContact++)
	{
		// vertex is from the v side, so the vertex order is always 3
		CFloatingType dis = accumulateVFContactForceAndHessian(iContact, 3, &contactResult, true, true, force, hessian);
		if (hessian.hasNaN())
		{
			std::cout << "hessian has NaN at vertex " << vId << std::endl;
			const auto fId = contactResult.contactPts[iContact].contactFaceId;
			const auto& mesh = getMesh(contactResult.contactPts[iContact].contactFaceSideMeshId);
			std::cout << "colliding with face: (" << mesh.facePosVId(fId, 0) << " " << mesh.facePosVId(fId, 1) << " " << mesh.facePosVId(fId, 2) 
				<< ") | distance: " << dis
				<< std::endl;
			std::exit(1);
		}
	}

	// neighbor faces collision detection
	for (size_t iNeiF = 0; iNeiF < pMesh->numNeiFaces(vId); iNeiF++)
	{
		IdType fId = pMesh->getVertexIthNeiFace(vId, iNeiF);

		CFloatingType faceMinDisToVertices = pClothContactDetector->getFaceMinDis(iMesh, fId);

		auto& faceContactInfos = pClothContactDetector->getFaceContactInfo(iMesh, fId);

		for (size_t iContact = 0; iContact < faceContactInfos.size(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			auto& faceContactInfo = faceContactInfos[iContact];
			ClothVFContactQueryResult& contactResultFV
				= getVFContactResult(faceContactInfo.meshIdVSide, faceContactInfo.vertexId);
			// getSimulatedMesh(faceContactInfo.meshIdVSide).clothClothVFContactQueyResults[faceContactInfo.vertexId];

			CFloatingType dis = accumulateVFContactForceAndHessian(faceContactInfo.contactId, pMesh->getVertexIthNeiFaceOrder(vId, iNeiF),
				&contactResultFV, true, true, force, hessian);
			//dMin = std::min(dis, dMin);
		}
	}

	// neighbor edges collision detection

	for (size_t iNeiE = 0; iNeiE < pMesh->numNeiEdges(vId); iNeiE++)
	{
		IdType eId = pMesh->getVertexIthNeiEdge(vId, iNeiE);
		ClothEEContactQueryResult& contactResultEE = getEEContactResult(iMesh, eId);

		for (size_t iContact = 0; iContact < contactResultEE.numContactPoints(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			CFloatingType dis = accumulateEEContactForceAndHessian(iContact, pMesh->getVertexIthNeiEdgeOrder(vId, iNeiE), &contactResultEE,
				true, true, force, hessian);
			// dMin = std::min(dis, dMin);
			if (hessian.hasNaN())
			{
				std::cout << "hessian has NaN at vertex " << vId << std::endl;
				const auto eId2 = contactResultEE.contactPts[iContact].contactEdgeId2;
				std::cout << "colliding between edge 1: (" << eId << "-(" << pMesh->getEdgeInfo(eId).eV1 << ", " << pMesh->getEdgeInfo(eId).eV2 << ")"
					<< std::endl;
				std::cout << " and between edge: " << eId2 << "-(" << pMesh->getEdgeInfo(eId2).eV1 << ", " << pMesh->getEdgeInfo(eId2).eV2 << std::endl << ")"
					<< std::endl;
				std::cout << "distance: " << dis << std::endl;;
 				const auto& mesh = getMesh(contactResultEE.contactPts[iContact].contactMeshId2);
				std::exit(1);
			}
		}
	}

// #define DEBUG_COLLISION_F_H

#ifdef DEBUG_COLLISION_F_H
	Vec3 col_force = force;
	Mat3 col_hessian = hessian;
#endif // DEBUG_COLLISION_F_H

	// accumulate the inertia and elasticity force and hessian
	pMesh->accumlateInertiaForceAndHessian(vId, force, hessian);
	pMesh->accumlateMaterialForceAndHessian(vId, force, hessian);
	if (pMesh->sewingVertices.size() > 0) {
		pMesh->accumlateSewingForceAndHessian(vId, force, hessian);
	}
	accumlateBoundaryForceAndHessian(pMesh, iMesh, vId, force, hessian);

	// accumulate the constraints' force and hessian
	accumulateConstraintForceAndHessian(iMesh, vId, force, hessian);

	// add external force
	force += pMesh->vertexExternalForces.col(vId);

	// Solve linear system
	Vec3 descentDirection;
	FloatingType stepSize = physicsParams().stepSize;
	bool solverSuccess = CuMatrix::solve3x3_psd_stable(hessian.data(), force.data(), descentDirection.data());
	if (descentDirection.hasNaN()) {
		std::cout << "descentDirection has NaN at vertex " << vId << std::endl;
#ifdef DEBUG_COLLISION_F_H
		std::cout
			<< "collision force: \n" << col_force << "\n"
			<< "collision hessian: \n" << col_hessian << "\n"
			<< "overall force: \n" << force << "\n"
			<< "overall hessian: \n" << hessian << "\n"
			<< std::endl;
		return;
#else
		std::exit(1);
#endif // DEBUG_COLLISION_F_H
	}
	if (!solverSuccess)
	{
		std::cout << "Solver failed at vertex " << vId << std::endl;
	}
	else
	{
		Vec3 positionNew = pMesh->vertex(vId) + stepSize * descentDirection;

		bool truncationHappened = applyConservativeBoundTruncation(iMesh, vId, positionNew);
		collisionDetectionRequired = truncationHappened;

		pMesh->vertexPositionNext(vId) = positionNew;
	}
}

void GAIA::VBDClothSimulationFramework::evaluateVertexForceAndHessian(IdType iMesh, IdType vId, Vec3& force, Mat3& hessian)
{
	VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);

	hessian = Mat3::Zero();
	force = Vec3::Zero();

	// vertex collision detection
	ClothVFContactQueryResult& contactResult = getVFContactResult(iMesh, vId);

	for (size_t iContact = 0; iContact < contactResult.numContactPoints(); iContact++)
	{
		// vertex is from the v side, so the vertex order is always 3
		CFloatingType dis = accumulateVFContactForceAndHessian(iContact, 3, &contactResult, true, true, force, hessian);
	}

	// neighbor faces collision detection
	for (size_t iNeiF = 0; iNeiF < pMesh->numNeiFaces(vId); iNeiF++)
	{
		IdType fId = pMesh->getVertexIthNeiFace(vId, iNeiF);

		CFloatingType faceMinDisToVertices = pClothContactDetector->getFaceMinDis(iMesh, fId);

		auto& faceContactInfos = pClothContactDetector->getFaceContactInfo(iMesh, fId);

		for (size_t iContact = 0; iContact < faceContactInfos.size(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			auto& faceContactInfo = faceContactInfos[iContact];
			ClothVFContactQueryResult& contactResultFV
				= getVFContactResult(faceContactInfo.meshIdVSide, faceContactInfo.vertexId);
			// getSimulatedMesh(faceContactInfo.meshIdVSide).clothClothVFContactQueyResults[faceContactInfo.vertexId];

			CFloatingType dis = accumulateVFContactForceAndHessian(faceContactInfo.contactId, pMesh->getVertexIthNeiFaceOrder(vId, iNeiF),
				&contactResultFV, true, true, force, hessian);
			//dMin = std::min(dis, dMin);
		}
	}

	// neighbor edges collision detection

	for (size_t iNeiE = 0; iNeiE < pMesh->numNeiEdges(vId); iNeiE++)
	{
		IdType eId = pMesh->getVertexIthNeiEdge(vId, iNeiE);
		ClothEEContactQueryResult& contactResultEE = getEEContactResult(iMesh, eId);

		for (size_t iContact = 0; iContact < contactResultEE.numContactPoints(); iContact++)
		{
			// vertex is from the f side, so the vertex order is pMesh->getVertexIthNeiFaceOrder(vId, iNeiF)
			CFloatingType dis = accumulateEEContactForceAndHessian(iContact, pMesh->getVertexIthNeiEdgeOrder(vId, iNeiE), &contactResultEE,
				true, true, force, hessian);
			// dMin = std::min(dis, dMin);
		}
	}

	// accumulate the inertia and elasticity force and hessian
	pMesh->accumlateInertiaForceAndHessian(vId, force, hessian);
	pMesh->accumlateMaterialForceAndHessian(vId, force, hessian);
	accumlateBoundaryForceAndHessian(pMesh, iMesh, vId, force, hessian);

	// accumulate the constraints' force and hessian
	accumulateConstraintForceAndHessian(iMesh, vId, force, hessian);

	// add external force
	force += pMesh->vertexExternalForces.col(vId);
}

bool GAIA::VBDClothSimulationFramework::applyConservativeBoundTruncation(IdType iMesh, IdType vId, Vec3& newPos)
{
	const VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
	Vec3 accumulatedDisplacement = newPos - pMesh->positionsAtPrevCollisionDetection.col(vId);
	FloatingType accumulatedDisplacementSize = accumulatedDisplacement.norm();

	CFloatingType convservativeBounds = pMesh->vertexConvervativeBounds[vId];

	if (accumulatedDisplacementSize > convservativeBounds)
		// moved out of the conservative bounds, truncation needed, and collision detection is required
	{
		FloatingType accumulatedDisplacementSizeNew = convservativeBounds;
		accumulatedDisplacement = accumulatedDisplacement * (accumulatedDisplacementSizeNew / accumulatedDisplacementSize);
		newPos = pMesh->positionsAtPrevCollisionDetection.col(vId) + accumulatedDisplacement;

		//std::cout << "Step truncation triggered for mesh " << iMesh << ", vertex " << vId << ", from " << accumulatedDisplacementSize	
		return true;
	}
	else
	{
		return false;

	}
}

CFloatingType GAIA::VBDClothSimulationFramework::accumulateVFContactForceAndHessian(IdType contactId, IdType contactVertexOrder,
	ClothVFContactQueryResult* pCollisionResult, bool updateContactInfo, bool apply_friction, Vec3& force, Mat3& hessian)
{
	VFContactPointInfo& contactPt = pCollisionResult->contactPts[contactId];

	FloatingType b;
	CFloatingType k = physicsParams().contactStiffness;

	if (updateContactInfo)
	{
		updateVFContactPointInfo(triMeshesAll, contactPt);
	}

	const TriMeshFEM* pMeshVSide = triMeshesAll[contactPt.contactVertexSideMeshId].get();
	const TriMeshFEM* pMeshFSide = triMeshesAll[contactPt.contactFaceSideMeshId].get();

	CVec3 x = pMeshVSide->vertex(contactPt.contactVertexId);
	CVec3 diff = x - contactPt.contactPoint; // from contact point to vertex
	CVec3 n = diff.normalized();
	FloatingType dis = diff.dot(n);

	if (dis < physicsParams().contactRadius + physicsParams().thickness)
	{
		CFloatingType bs[4] = { -contactPt.barycentrics[0], -contactPt.barycentrics[1], -contactPt.barycentrics[2], 1.f };
		CFloatingType b = bs[contactVertexOrder];

		// simplified point-plane
		// Penalty force
		FloatingType dEdD, d2E_dDdD;
		computeContactRepulsiveForce(dis, physicsParams().thickness, dEdD, d2E_dDdD);
		CFloatingType lambda = -dEdD;

		force += b * lambda * n;
		hessian += d2E_dDdD * b * b * n * n.transpose();

		// Friction
		if (apply_friction) {
			IdType t0 = pMeshFSide->facePosVId(contactPt.contactFaceId, 0);
			IdType t1 = pMeshFSide->facePosVId(contactPt.contactFaceId, 1);
			IdType t2 = pMeshFSide->facePosVId(contactPt.contactFaceId, 2);
			CVec3 dx3 = x - pMeshVSide->vertexPrevPos(contactPt.contactVertexId);
			CVec3& bary = contactPt.barycentrics;
			CVec3 contactPtPrev = bary.x() * pMeshFSide->vertexPrevPos(t0)
				+ bary.y() * pMeshFSide->vertexPrevPos(t1)
				+ bary.z() * pMeshFSide->vertexPrevPos(t2);
			//Vec3 dx0 = pMeshFSide->vertex(t0) - pMeshFSide->vertexPrevPos(t0);
			//Vec3 dx1 = pMeshFSide->vertex(t1) - pMeshFSide->vertexPrevPos(t1);
			//Vec3 dx2 = pMeshFSide->vertex(t2) - pMeshFSide->vertexPrevPos(t2);
			//Vec3 dx = dx3 - (bary.x() * dx0 + bary.y() * dx1 + bary.z() * dx2);
			Vec3 dx = dx3 - (contactPt.contactPoint - contactPtPrev);
			Vec3 e0 = (pMeshFSide->vertex(t1) - pMeshFSide->vertex(t0)).normalized();
			Mat3x2 T = Mat3x2::Zero();
			T.col(0) = e0;
			T.col(1) = (e0.cross(n)).normalized();
			Vec2 u = T.transpose() * dx;

			// average of the two friction coefficients
			CFloatingType mu = (pMeshFSide->pObjectParams->frictionDynamic
				+ pMeshVSide->pObjectParams->frictionDynamic) * 0.5;
			CFloatingType epsV = (pMeshFSide->pObjectParams->frictionEpsV
				+ pMeshVSide->pObjectParams->frictionEpsV) * 0.5;
			CFloatingType dt = physicsParams().dt;
			CFloatingType epsU = epsV * dt;
			Vec3 frictionForce;
			Mat3 frictionHessian;
			computeFriction(mu, lambda, T, u, epsU, frictionForce, frictionHessian);

			force += b * frictionForce;
			hessian += b * b * frictionHessian;
		}
	}

	return dis;
}

CFloatingType GAIA::VBDClothSimulationFramework::accumulateEEContactForceAndHessian(IdType contactId, IdType contactVertexOrder,
	ClothEEContactQueryResult* pCollisionResult, bool updateContactInfo, bool apply_friction, Vec3& force, Mat3& hessian)
{
	EEContactPointInfo& eeContactPt = pCollisionResult->contactPts[contactId];
	if (updateContactInfo)
	{
		updateEEContactPointInfo(triMeshesAll, eeContactPt);
	}

	const Vec3 x1 = eeContactPt.c1;
	const Vec3 x2 = eeContactPt.c2;

	CVec3 diff = x1 - x2; // from contact point 2 to contact point 1
	CFloatingType dis = diff.norm();

	const TriMeshFEM* pMesh1 = triMeshesAll[eeContactPt.contactMeshId1].get();
	const TriMeshFEM* pMesh2 = triMeshesAll[eeContactPt.contactMeshId2].get();
	const EdgeInfo& edgeInfo1 = pMesh1->pTopology->edgeInfos[eeContactPt.contactEdgeId1];
	const EdgeInfo& edgeInfo2 = pMesh2->pTopology->edgeInfos[eeContactPt.contactEdgeId2];

	CVec3 v1 = pMesh1->vertex(edgeInfo1.eV2) - pMesh1->vertex(edgeInfo1.eV1);
	CVec3 v2 = pMesh2->vertex(edgeInfo2.eV2) - pMesh2->vertex(edgeInfo2.eV1);

	CFloatingType parallelEps = v1.cross(v2).squaredNorm() / (v1.norm() * v2.norm());

	if (dis < physicsParams().contactRadius + physicsParams().thickness
		// && dis > 1e-7f 
		&& parallelEps > CMP_EPSILON
		)
	{
		CVec3 n = diff / dis;

		FloatingType bs[4] = { 1.f - eeContactPt.mu1, eeContactPt.mu1, -1.f + eeContactPt.mu2, -eeContactPt.mu2 };
		FloatingType b = bs[contactVertexOrder];
		CFloatingType k = physicsParams().contactStiffness;
		CFloatingType penetrationDepth = physicsParams().contactRadius - dis;

		FloatingType dEdD, d2E_dDdD;
		computeContactRepulsiveForce(dis, physicsParams().thickness, dEdD, d2E_dDdD);
		CFloatingType lambda = -dEdD;
		force += b * lambda * n;
		hessian += d2E_dDdD * b * b * n * n.transpose();

		// Friction
		if (apply_friction) {
			CVec3 c1Prev = (1.f - eeContactPt.mu1) * pMesh1->vertexPrevPos(edgeInfo1.eV1)
				+ eeContactPt.mu1 * pMesh1->vertexPrevPos(edgeInfo1.eV2);
			CVec3 c2Prev = (1.f - eeContactPt.mu2) * pMesh2->vertexPrevPos(edgeInfo2.eV1)
				+ eeContactPt.mu2 * pMesh2->vertexPrevPos(edgeInfo2.eV2);
			//Vec3 dx0 = pMeshFSide->vertex(t0) - pMeshFSide->vertexPrevPos(t0);
			//Vec3 dx1 = pMeshFSide->vertex(t1) - pMeshFSide->vertexPrevPos(t1);
			//Vec3 dx2 = pMeshFSide->vertex(t2) - pMeshFSide->vertexPrevPos(t2);
			//Vec3 dx = dx3 - (bary.x() * dx0 + bary.y() * dx1 + bary.z() * dx2);
			Vec3 dx = (eeContactPt.c1 - c1Prev) - (eeContactPt.c2 - c2Prev);
			Mat3x2 T = Mat3x2::Zero();
			T.col(0) = v1.normalized();
			T.col(1) = (T.col(0).cross(n)).normalized();
			Vec2 u = T.transpose() * dx;

			// average of the two friction coefficients
			CFloatingType mu = (pMesh1->pObjectParams->frictionDynamic
				+ pMesh2->pObjectParams->frictionDynamic) * 0.5;
			CFloatingType epsV = (pMesh1->pObjectParams->frictionEpsV
				+ pMesh2->pObjectParams->frictionEpsV) * 0.5;
			CFloatingType dt = physicsParams().dt;
			CFloatingType epsU = epsV * dt;
			Vec3 frictionForce;
			Mat3 frictionHessian;
			computeFriction(mu, lambda, T, u, epsU, frictionForce, frictionHessian);

			force += b * frictionForce;
			hessian += b * b * frictionHessian;
		}
	}
	return dis;

}

void GAIA::VBDClothSimulationFramework::accumlateBoundaryForceAndHessian(VBDBaseTriMesh* pMesh, IdType iMesh, IdType vId, Vec3& force, Mat3& hessian)
{
	CFloatingType boundaryContactStiffness = physicsParams().boundaryContactStiffness;
	CFloatingType boundaryFrictionDynamic = physicsParams().boundaryFrictionDynamic;
	CFloatingType boundaryFrictionEpsValue = physicsParams().boundaryFrictionEpsValue;

	for (size_t iDim = 0; iDim < 3; iDim++)
	{
		Vec3 contactNormal = Vec3::Zero();
		CFloatingType dt = physicsParams().dt;

		CFloatingType lowerBound = physicsParams().worldBounds(iDim, 0);
		CFloatingType upperBound = physicsParams().worldBounds(iDim, 1);

		if (pMesh->vertex(vId)[iDim] < lowerBound)
		{
			// Penalty Force
			CFloatingType penetrationDepth = lowerBound - pMesh->vertex(vId)[iDim];
			force(iDim) += penetrationDepth * boundaryContactStiffness;
			hessian(iDim, iDim) += boundaryContactStiffness;

			// Friction
			Vec3 dx = pMesh->vertex(vId) - pMesh->vertexPrevPos(vId);
			Mat3x2 T = Mat3x2::Zero();
			T((iDim + 1) % 3, 0) = 1;
			T((iDim + 2) % 3, 1) = 1;
			Vec2 u = T.transpose() * dx;
			CFloatingType lambda = penetrationDepth * boundaryContactStiffness;
			CFloatingType mu = boundaryFrictionDynamic;
			CFloatingType epsV = boundaryFrictionEpsValue;
			CFloatingType epsU = epsV * dt;
			Vec3 frictionForce;
			Mat3 frictionForceHessian;
			computeFriction(mu, lambda, T, u, epsU, frictionForce, frictionForceHessian);

			force += frictionForce;
			hessian += frictionForceHessian;

		}
		else if (pMesh->vertex(vId)[iDim] > upperBound)
		{
			CFloatingType penetrationDepth = pMesh->vertex(vId)[iDim] - upperBound;
			force(iDim) -= penetrationDepth * physicsParams().boundaryContactStiffness;
			hessian(iDim, iDim) += physicsParams().boundaryContactStiffness;

			//// Friction
			Vec3 dx = pMesh->vertex(vId) - pMesh->vertexPrevPos(vId);
			Mat3x2 T = Mat3x2::Zero();
			T((iDim + 1) % 3, 0) = 1;
			T((iDim + 2) % 3, 1) = 1;
			Vec2 u = T.transpose() * dx;
			CFloatingType lambda = penetrationDepth * boundaryContactStiffness;
			CFloatingType mu = boundaryFrictionDynamic;
			CFloatingType epsV = boundaryFrictionEpsValue;
			CFloatingType epsU = epsV * dt;
			Vec3 frictionForce;
			Mat3 frictionForceHessian;
			computeFriction(mu, lambda, T, u, epsU, frictionForce, frictionForceHessian);
			force += frictionForce;
			hessian += frictionForceHessian;
		}
	}
}

void GAIA::VBDClothSimulationFramework::accumulateConstraintForceAndHessian(IdType iMesh, IdType vId, Vec3& force, Mat3& hessian)
{
	VertexConstraints& constraints = getSimulatedMesh(iMesh).vertexConstraints[vId];
	for (IdType iConstraint = 0; iConstraint < constraints.attachmentConstraints.size(); iConstraint++)
	{
		VertexFaceAttachmentConstraint attachmentConstraint = constraints.attachmentConstraints[iConstraint];

		Vec3 vertPos = getSimulatedMesh(attachmentConstraint.attachedMeshVSdie).vertex(vId);

		IdType attachedMeshIdFaceSide = attachmentConstraint.attachedMeshFaceSide;

		CFloatingType bs[4] = { -attachmentConstraint.attachedPos[0], -attachmentConstraint.attachedPos[1],
			-attachmentConstraint.attachedPos[2], 1.f };
		Vec3 attachedPos;
		if (attachedMeshIdFaceSide != -1)
			// attach to a face
		{
			TriMeshFEM* attachedMesh = &getMesh(attachedMeshIdFaceSide);
			attachedPos = attachmentConstraint.attachedPos[0] * attachedMesh->vertex(attachmentConstraint.attachedFaceVIds[0])
				+ attachmentConstraint.attachedPos[1] * attachedMesh->vertex(attachmentConstraint.attachedFaceVIds[1])
				+ attachmentConstraint.attachedPos[2] * attachedMesh->vertex(attachmentConstraint.attachedFaceVIds[2]);
		}
		else
			// attach to a fixed point
		{
			attachedPos = attachmentConstraint.attachedPos;
		}

		Vec3 diff = vertPos - attachedPos;
		CFloatingType b = bs[attachmentConstraint.vertexOrder];
		force -= attachmentConstraint.constraintStiffness * b * diff;
		hessian += attachmentConstraint.constraintStiffness * b * b * Mat3::Identity();
	}
}

void GAIA::VBDClothSimulationFramework::computeFriction(CFloatingType mu, CFloatingType lambda, const Mat3x2& T,
	const Vec2& u, CFloatingType epsU, Vec3& force, Mat3& hessian)
{
	// Friction
	CFloatingType uNorm = u.norm();
	if (uNorm > 0)
	{
		// IPC friction 
		// https://github.com/ipc-sim/ipc-toolkit/blob/main/src/ipc/friction/smooth_friction_mollifier.cpp
		FloatingType f1_SF_over_x;
		FloatingType df1_x_minus_f1_over_x3;
		if (uNorm > epsU)
		{
			f1_SF_over_x = 1 / uNorm;
			df1_x_minus_f1_over_x3 = -1 / (uNorm * uNorm * uNorm);
		}
		else
		{
			f1_SF_over_x = (-uNorm / epsU + 2) / epsU;
			df1_x_minus_f1_over_x3 = -1 / (uNorm * epsU * epsU);
		}
		force = -mu * lambda * T * f1_SF_over_x * u;
		// hessian = mu * lambda * T * (df1_x_minus_f1_over_x3 * u * u.transpose() + f1_SF_over_x * Mat2::Identity()) * T.transpose();
		// IPC made a mistake in differetiating through n, which makes it incredibly instable :)
		hessian = mu * lambda * T * (f1_SF_over_x * Mat2::Identity()) * T.transpose();
	}
	else
	{
		force.setZero();
		hessian.setZero();
	}
}

void GAIA::VBDClothSimulationFramework::computeContactRepulsiveForce(CFloatingType dis_, CFloatingType thickness, FloatingType& dEdD, FloatingType& d2E_dDdD)
{
	// switch (physicsParams().contactRadius)
	FloatingType disAdjusted = dis_ - thickness;
	CFloatingType contactRadius = physicsParams().contactRadius;
	CFloatingType penetrationDepth = contactRadius - disAdjusted;
	CFloatingType k = physicsParams().contactStiffness;


	switch (physicsParams().contactEnergyTypes)
	{
	case 0:
		dEdD = -k * penetrationDepth;
		d2E_dDdD = k;
		break;

	case 1:
		// 1 : quadratic - logrithmic - 2 - stages
	{
		CFloatingType tau = contactRadius * 0.5f; // to make energy C2 continuous
		if (disAdjusted < tau && disAdjusted > CMP_EPSILON)
		{
			//CFloatingType k2 = tau * k * (tau - physicsParams().contactRadius);
			CFloatingType k2 = 0.5f * SQR(tau) * k;
			dEdD = -k2 / disAdjusted;
			d2E_dDdD = k2 / SQR(disAdjusted);
		}
		else {
			dEdD = -k * penetrationDepth;
			d2E_dDdD = k;
		}

	}
	break;
	case 2:
		dEdD = (- k * SQR(disAdjusted - contactRadius) / disAdjusted - 2 * k * log(disAdjusted / contactRadius) * (disAdjusted - contactRadius));
		d2E_dDdD = k * ((SQR(disAdjusted - contactRadius) / SQR(disAdjusted)) - (4 * disAdjusted - 4 * contactRadius) / disAdjusted - 2 * log(disAdjusted / contactRadius));

		if (dEdD == NAN || d2E_dDdD == NAN)
		{
			printf("!");
		}
	break;
	default:
		break;
	}
}

void GAIA::VBDClothSimulationFramework::computeContactEnergy(CFloatingType dis, CFloatingType thickness, FloatingType& E)
{
	// switch (physicsParams().contactRadius)
	FloatingType disAdjusted = dis - thickness;
	CFloatingType contactRadius = physicsParams().contactRadius;
	CFloatingType penetrationDepth = contactRadius - disAdjusted;
	CFloatingType k = physicsParams().contactStiffness;


	switch (physicsParams().contactEnergyTypes)
	{
	case 0:
		E = 0.5f * k * SQR(penetrationDepth);
		break;

	case 1:
		// 1 : quadratic - logrithmic - 2 - stages
	{
		CFloatingType tau = contactRadius * 0.5f; // to make energy C2 continuous
		if (disAdjusted < tau && disAdjusted > CMP_EPSILON)
		{
			//CFloatingType k2 = tau * k * (tau - physicsParams().contactRadius);
			CFloatingType k2 = 0.5f * SQR(tau) * k;
			CFloatingType b = k2 / contactRadius + k2 * log(tau);
			E = -log(disAdjusted) * k2 + b;
		}
		else {
			E = 0.5f * k * SQR(penetrationDepth);
		}

	}
	break;
	default:
		E = 0;
		break;
	}
}

//void GAIA::VBDClothSimulationFramework::simulate()
//{
//	if (basePhysicsParams->checkAndUpdateWorldBounds)
//	{
//		updateWorldBox();
//	}
//	setUpOutputFolders(outputFolder);
//	std::cout
//		<< "----------------------------------------------------\n"
//		<< "Output folder is: " << outputFolder << std::endl;
//
//	writeOutputs(outputFolder, frameId);
//	std::cout
//		<< "----------------------------------------------------\n"
//		<< "Starting Sims\n"
//		<< "----------------------------------------------------"
//		<< std::endl;
//
//	baseTimeStatistics->setToZero();
//
//	while (frameId < basePhysicsParams->numFrames) {
//		TICK(timeCsmpFrame);
//
//		debugOperation(DEBUG_LVL_INFO, [&]() {
//			std::cout
//				<< "----------------------------------------------------\n"
//				<< "Frame " << frameId + 1 << " begin.\n"
//				;
//			});
//
//		// loadTargetMesh(frameId);
//
//		/*std::vector<int> breach = { 23, 73, 124, 175, 226, 277, 326, 375, 424, 473,
//			524, 575, 626, 677, 728, 779, 830, 879, 928, 977, 1026, 1075, 1124, 1175, 
//			1226, 1277, 1328, 1379, 1430, 1479, 1528, 1577, 1626, 1675, 1724, 1773, 1822,
//			1871, 1920, 1971, 2022, 2073, 2124, 2175, 2226, 2277, 2327, 2377, 2428, 2478 };
//		int startFrame = 240;
//		int tearGap = 3;
//		if (frameId >= startFrame)
//		{
//			if (frameId % tearGap == 0) {
//				int i = (frameId-startFrame) / tearGap;
//
//				if (i<breach.size()-1)
//				{
//					tearMesh(0, breach[i], breach[i+1]
//						);
//				}
//
//			}
//		}*/
//
//		if (physicsParams().useParallel)
//		{
//			runStep();
//		}
//		else {
//			runStepSequential();
//		}
//
//		TICK(timeCsmpSaveOutputs);
//		writeOutputs(outputFolder, frameId + 1);
//		TOCK_STRUCT((*baseTimeStatistics), timeCsmpSaveOutputs);
//
//		TOCK_STRUCT((*baseTimeStatistics), timeCsmpFrame);
//
//		debugPrint(DEBUG_LVL_INFO, baseTimeStatistics->getString());
//		debugOperation(DEBUG_LVL_INFO, [&]() {
//			std::cout
//				<< "Frame " << frameId + 1 << " completed, Time consumption: " << baseTimeStatistics->timeCsmpFrame << "\n"
//				<< "----------------------------------------------------\n";
//			});
//
//		++frameId;
//		baseTimeStatistics->setToZero();
//	}
//}

void GAIA::VBDClothSimulationFramework::tearMesh(IdType iMesh, IdType v1, IdType v2)
{
	auto& mesh = getSimulatedMesh(iMesh);
	int ret = mesh.tearMesh(v1, v2);
	// update coloring
	// ret = 0: no duplicate vertices
	// ret = 1: duplicate v1
	// ret = 2: duplicate v2
	// ret = 3: duplicate v1 and v2
	if (ret == 1) {
		IdType v3 = mesh.numVertices() - 1;
		auto color = mesh.globalColors[v1];
		vertexParallelGroups[color].push_back(iMesh);
		vertexParallelGroups[color].push_back(v3);
	}
	else if (ret == 2) {
		IdType v4 = mesh.numVertices() - 1;
		auto color = mesh.globalColors[v2];
		vertexParallelGroups[color].push_back(iMesh);
		vertexParallelGroups[color].push_back(v4);
	}
	else if (ret == 3) {
		IdType v3 = mesh.numVertices() - 2;
		IdType v4 = mesh.numVertices() - 1;
		auto color = mesh.globalColors[v1];
		vertexParallelGroups[color].push_back(iMesh);
		vertexParallelGroups[color].push_back(v3);
		color = mesh.globalColors[v2];
		vertexParallelGroups[color].push_back(iMesh);
		vertexParallelGroups[color].push_back(v4);
	}
}

void GAIA::VBDClothSimulationFramework::applyInitialGuess()
{
	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		mesh.applyInitialStep();

		mesh.evaluateExternalForce();
	}
}

FloatingType GAIA::VBDClothSimulationFramework::evaluateEnergy(ConvergenceStats& energyStats)
{
	energyStats.energy_allClothes = 0.f;
	energyStats.avgGradNorm_allClothes = 0.f;
	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);

		FloatingType me, me_inertia, me_elastic, me_elastic_StVK, me_elastic_bending;
		me = mesh.evaluatedMeritEnergy(me_inertia, me_elastic_StVK, me_elastic_bending);
		me_elastic = me_elastic_StVK + me_elastic_bending;

		if (!physicsParams().withInertia)
		{
			me = me_elastic;
			me_inertia = 0.f;
		}

		energyStats.energy_allClothes += me;
		energyStats.energy[iCloth] = me;
		energyStats.energy_inertia[iCloth] = me_inertia;
		energyStats.energy_elastic[iCloth] = me_elastic;
		energyStats.energy_elastic_StVK[iCloth] = me_elastic_StVK;
		energyStats.energy_elastic_bending[iCloth] = me_elastic_bending;

		/*if (physicsParams().fit)
		{
			FloatingType energyReg = mesh.evaluatedRegistrationEnergy();
			energyStats.energy_registration[iCloth] = energyReg;
			energyStats.energy_allClothes += energyReg;
		}
		else {
			energyStats.energy_registration[iCloth] = 0.f;
		}

		FloatingType avgGradNorm = mesh.evaluatedAvgGradNorm();
		energyStats.avgGradNorm[iCloth] = mesh.evaluatedAvgGradNorm();
		energyStats.avgGradNorm_allClothes += avgGradNorm;*/
	}

	energyStats.avgGradNorm_allClothes = energyStats.avgGradNorm_allClothes / numClothes();

	energyStats.energy_allClothes += evaluateCollisionEnergy(energyStats);

	return energyStats.energy_allClothes;
}

FloatingType GAIA::VBDClothSimulationFramework::evaluateCollisionEnergy(ConvergenceStats& energyStats)
{
	FloatingType collisionEnergyAll = 0.f;
	const FloatingType collisionStiffness = physicsParams().contactStiffness;
	const FloatingType r = pClothContactDetectorParameters->maxQueryDis;
	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		//std::atomic<FloatingType> collisionEnergy = 0.f;
		FloatingType collisionEnergy = 0.f;

		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		switch (physicsParams().contactEnergyTypes)
		{
			// 0: quadratic
			for (size_t iVert = 0; iVert < mesh.numVertices(); iVert++)
			{
				const ClothVFContactQueryResult& contactPointResult = getVFContactResult(iCloth, iVert);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						FloatingType d = (contactPointResult.contactPts[iContact].contactPoint - mesh.vertex(iVert)).norm();
						if (d < r)
						{
							collisionEnergy += 0.5f * collisionStiffness * (r - d) * (r - d);
						}
					}
				}
			}
		case 1:
			// 1: logrithmic
			//cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) 
			for (size_t iVert = 0; iVert < mesh.numVertices(); iVert++)
			{
				const ClothVFContactQueryResult& contactPointResult = getVFContactResult(iCloth, iVert);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						FloatingType d = (contactPointResult.contactPts[iContact].contactPoint - mesh.vertex(iVert)).norm();
						if (d < r)
						{
							collisionEnergy += collisionStiffness * log(d) - log(r);
						}
					}
				}
			}
			//);
		default:
			break;
		}

		switch (physicsParams().contactEnergyTypes)
		{
		case 0:
			// 0: quadratic
			for (size_t iE = 0; iE < mesh.numEdges(); iE++)
			{
				const ClothEEContactQueryResult& contactPointResult = getEEContactResult(iCloth, iE);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						const auto& contactPt = contactPointResult.contactPts[iContact];
						const Vec3 C1 = mesh.edgeLerp(iE, contactPt.mu1);
						// the contacting mesh & edge from the other side
						const VBDBaseTriMesh& meshContacting = getSimulatedMesh(contactPt.contactMeshId2);
						const Vec3 C2 = meshContacting.edgeLerp(contactPt.contactEdgeId2, contactPt.mu2);
						const EdgeInfo& edgeInfo1 = mesh.getEdgeInfo(iE);
						const FloatingType d = (C2 - C1).norm();

						if (d < r)
						{
							collisionEnergy += 0.5f * collisionStiffness * (r - d) * (r - d);

						}
					}
				}
			}
			break;
		case 1:
			// 1: logrithmic
			for (size_t iE = 0; iE < mesh.numEdges(); iE++)
			{
				const ClothEEContactQueryResult& contactPointResult = getEEContactResult(iCloth, iE);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						const auto& contactPt = contactPointResult.contactPts[iContact];
						const Vec3 C1 = mesh.edgeLerp(iE, contactPt.mu1);
						const VBDBaseTriMesh& meshContacting = getSimulatedMesh(contactPt.contactMeshId2);
						const Vec3 C2 = meshContacting.edgeLerp(contactPt.contactEdgeId2, contactPt.mu2);

						const FloatingType d = (C2 - C1).norm();

						if (d < r)
						{
							collisionEnergy += collisionStiffness * log(d) - log(r);
						}
					}
				}
			}
		default:
			break;
		}

		energyStats.energy_collision[iCloth] = collisionEnergy;
		collisionEnergyAll += collisionEnergy;
	}

	return collisionEnergyAll;
}

void GAIA::VBDClothSimulationFramework::clearGradient()
{
	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		mesh.clearGradient();
	}
}

void GAIA::VBDClothSimulationFramework::evaluateConvergence()
{
	FloatingType avgGradNorm = 0.f;
	FloatingType maxGradNorm = 0.f;

	if (substep + 1 > runtimeStatistics().avgForceResiduals.size() && iIter == 0)
	{
		runtimeStatistics().avgForceResiduals.emplace_back();
		runtimeStatistics().maxForceResiduals.emplace_back();
	}

	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iMesh);
		mesh.gradient.setZero();
		cpu_parallel_for(0, mesh.numVertices(), [&](int vId) {
			//for (size_t iV = 0; iV < numVertices; iV++) {
			auto pMesh = &getSimulatedMesh(iMesh);
			Vec3 force = Vec3::Zero();
			Mat3 hessian = Mat3::Zero();
			if (!pMesh->fixedMask[vId])
			{
				evaluateVertexForceAndHessian(iMesh, vId, force, hessian);
			}
			mesh.gradient.col(vId) = force;
		});

		for (size_t iV = 0; iV < mesh.numVertices(); iV++)
		{
			avgGradNorm += mesh.gradient.col(iV).norm();
			maxGradNorm = std::max(maxGradNorm, mesh.gradient.col(iV).norm());
		}
	}
	avgGradNorm = avgGradNorm / numAllVertices;

	runtimeStatistics().avgForceResiduals.back().push_back(avgGradNorm);
	runtimeStatistics().maxForceResiduals.back().push_back(maxGradNorm);
}

//void GAIA::VBDClothSimulationFramework::preconditionGradients()
//{
//
//	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
//	{
//		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//		if (physicsParams().usePreconditioner)
//		{
//			cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
//				// assemble the diagonal block of the Hessian matrix
//				// first step, the internal force
//				Mat3 h = mesh.vertexInternalForces_HessianDiagonalBlocks[iVert];
//				Vec3 gradPreconditioned;
//				const Vec3 grad = mesh.gradient.col(iVert);
//				// then add the collision force
//
//				// add the registration force
//
//				CuMatrix::solve3x3_psd_stable(h.data(), grad.data(), gradPreconditioned.data());
//				// gradPreconditioned = h.colPivHouseholderQr().solve(grad);
//
//				if (gradPreconditioned.dot(grad) > CMP_EPSILON)
//				{
//					mesh.gradientPreconditioned.col(iVert) = gradPreconditioned;
//				}
//				else
//				{
//					if (grad.norm() != 0)
//					{
//						//std::cout << "Warning: negative gradient after preconditioning, revert to original gradient." << std::endl;
//					}
//					mesh.gradientPreconditioned.col(iVert) = grad;
//				}
//				});
//		}
//		else
//		{
//			mesh.gradientPreconditioned = mesh.gradient;
//		}
//	}
//}

//void GAIA::VBDClothSimulationFramework::accumulateMaterialGradient()
//{
//	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
//	{
//		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//		mesh.accumulateMaterialGradient(physicsParams().withInertia);
//
//	}
//}

void GAIA::VBDClothSimulationFramework::accumulateRegistrationGradient()
{
	if (!physicsParams().fit)
	{
		return;
	}
	if (iIter == 0 || !(iIter % physicsParams().closestPointQuerySteps))
	{
		queryClosestTargetPoint();
	}

	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
			/*const TriMeshClosestPointQueryResult& targetPointResult = mesh.targetPointsQueyResults[iVert];
			if (targetPointResult.found)
			{
				mesh.gradient.col(iVert) +=
					physicsParams().registrationEnergyWeight * (mesh.vertex(iVert) - targetPointResult.contactPts.back().closestPt);
			}*/
			});
	}
}

void GAIA::VBDClothSimulationFramework::accumulateCollisionGradient()
{
	if (contactInfoUpdated)
	{
		contactInfoUpdated = false;
	}
	else {
		recomputeContacts();
	}

	const FloatingType r = pClothContactDetectorParameters->maxQueryDis;
	for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);

		switch (physicsParams().contactEnergyTypes)
		{
		case 0:
			// 0: quadratic
			cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
				const ClothVFContactQueryResult& contactPointResult = getVFContactResult(iCloth, iVert);
				if (contactPointResult.found && !mesh.fixedMask(iVert))
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						const Vec3 diff = contactPointResult.contactPts[iContact].contactPoint - mesh.vertex(iVert);
						const FloatingType d = contactPointResult.contactPts[iContact].d;
						if (d < r)
						{
							mesh.gradient.col(iVert) +=
								physicsParams().contactStiffness * diff * (r - d) / d;
						}
					}
				}
				});
			break;
		case 1:
			// 1: logrithmic
			cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
				const ClothVFContactQueryResult& contactPointResult = getVFContactResult(iCloth, iVert);
				if (contactPointResult.found && !mesh.fixedMask(iVert))
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						FloatingType d = contactPointResult.contactPts[iContact].d;
						if (d < r)
						{
							mesh.gradient.col(iVert) +=
								physicsParams().contactStiffness * (contactPointResult.contactPts[iContact].contactPoint - mesh.vertex(iVert)) / (d * d);
						}
					}
				}
				});
			break;
		default:
			break;
		}

		// EE contact

		switch (physicsParams().contactEnergyTypes)
		{
		case 0:
			// 0: quadratic
			for (size_t iE = 0; iE < mesh.numEdges(); iE++)
			{
				const ClothEEContactQueryResult& contactPointResult = getEEContactResult(iCloth, iE);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						const auto& contactPt = contactPointResult.contactPts[iContact];
						const Vec3& C1 = contactPt.c1;
						const Vec3& C2 = contactPt.c2;
						const EdgeInfo& edgeInfo1 = mesh.getEdgeInfo(iE);
						const FloatingType d = contactPt.d;

						const Vec3 grad = physicsParams().contactStiffness * (C2 - C1) * (r - d) / d;

						FloatingType t1 = (1.f - contactPt.mu1);

						const IdType eV1 = edgeInfo1.eV1;
						const IdType eV2 = edgeInfo1.eV2;

						if (d < r)
						{
							if (!mesh.fixedMask(eV1))
							{
								mesh.gradient.col(eV1) += t1 * grad; // *mesh.vertexInvMass(eV1);
							}
							if (!mesh.fixedMask(eV2))
							{
								mesh.gradient.col(eV2) += contactPt.mu1 * grad;// *mesh.vertexInvMass(eV2);
							}
						}
					}
				}
			}
			break;
		case 1:
			// 1: logrithmic
			for (size_t iE = 0; iE < mesh.numEdges(); iE++)
			{
				const ClothEEContactQueryResult& contactPointResult = getEEContactResult(iCloth, iE);
				if (contactPointResult.found)
				{
					for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
					{
						const auto& contactPt = contactPointResult.contactPts[iContact];
						const Vec3& C1 = contactPt.c1;
						const Vec3& C2 = contactPt.c2;
						const EdgeInfo& edgeInfo1 = mesh.getEdgeInfo(iE);
						const FloatingType d = contactPt.d;

						const Vec3 grad = physicsParams().contactStiffness * (C2 - C1) / (d * d);

						FloatingType t1 = (1.f - contactPt.mu1);

						const IdType eV1 = edgeInfo1.eV1;
						const IdType eV2 = edgeInfo1.eV2;

						if (d < r)
						{
							if (!mesh.fixedMask(eV1))
							{
								mesh.gradient.col(eV1) += t1 * grad; // *mesh.vertexInvMass(eV1);
							}
							if (!mesh.fixedMask(eV2))
							{
								mesh.gradient.col(eV2) += contactPt.mu1 * grad; // *mesh.vertexInvMass(eV2);
							}
						}
					}
				}
			}
		default:
			break;
		}

	}
}

void GAIA::VBDClothSimulationFramework::recomputeContacts()
{
	for (int iCloth = 0; iCloth < numClothes(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
			// recompute VF contacts
			ClothVFContactQueryResult& contactPointResult = getVFContactResult(iCloth, iVert);

			for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
			{
				auto& contactPt = contactPointResult.contactPts[iContact];

				const VBDBaseTriMesh& meshFaceSide = getSimulatedMesh(contactPt.contactFaceSideMeshId);
				const IdType* face = mesh.facePos.col(contactPt.contactFaceId).data();

				const embree::Vec3fa a = embree::Vec3fa::loadu(meshFaceSide.vertex(face[0]).data());
				const embree::Vec3fa b = embree::Vec3fa::loadu(meshFaceSide.vertex(face[1]).data());
				const embree::Vec3fa c = embree::Vec3fa::loadu(meshFaceSide.vertex(face[2]).data());

				ClosestPointOnTriangleType pointType;
				embree::Vec3fa closestPtBarycentrics;
				const embree::Vec3fa queryPt = embree::Vec3fa::loadu(mesh.vertex(iVert).data());

				const embree::Vec3fa closestP = GAIA::closestPointTriangle(queryPt, a, b, c, closestPtBarycentrics, pointType);
				float d = embree::distance(queryPt, closestP);
				contactPt.d = d;
				contactPt.contactPoint << closestP.x, closestP.y, closestP.z;
				contactPt.closestPtType = pointType;
			}
			});

		cpu_parallel_for(0, mesh.numEdges(), [&](int iE) {
			// recompute EE contacts
			ClothEEContactQueryResult& contactPointResult = getEEContactResult(iCloth, iE);

			const EdgeInfo& edgeInfoQuery = mesh.pTopology->edgeInfos[iE];

			for (size_t iContact = 0; iContact < contactPointResult.contactPts.size(); iContact++)
			{
				auto& contactPt = contactPointResult.contactPts[iContact];
				const VBDBaseTriMesh& meshContacting = getSimulatedMesh(contactPt.contactMeshId2);
				const EdgeInfo& edgeInfoTarget = meshContacting.pTopology->edgeInfos[contactPt.contactEdgeId2];

				const embree::Vec3fa p1 = embree::Vec3fa::loadu(mesh.vertex(edgeInfoQuery.eV1).data());
				const embree::Vec3fa p2 = embree::Vec3fa::loadu(mesh.vertex(edgeInfoQuery.eV2).data());

				const embree::Vec3fa q1 = embree::Vec3fa::loadu(meshContacting.vertex(edgeInfoTarget.eV1).data());
				const embree::Vec3fa q2 = embree::Vec3fa::loadu(meshContacting.vertex(edgeInfoTarget.eV2).data());
				embree::Vec3fa c1, c2;
				FloatingType mua, mub;
				get_closest_points_between_segments(p1, p2, q1, q2, c1, c2, mua, mub);

				FloatingType d = embree::distance(c1, c2);

				contactPt.mu1 = mua;
				contactPt.mu2 = mub;
				contactPt.d = d;
				contactPt.c1 << c1.x, c1.y, c1.z;
				contactPt.c2 << c2.x, c2.y, c2.z;
			}

			});
	}

	contactInfoUpdated = true;
}

void GAIA::VBDClothSimulationFramework::queryClosestTargetPoint()
{
	for (IdType iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		//cpu_parallel_for(0, mesh.numVertices(), [&](int iVert) {
		//	pTargetMatcher->closestPointQuery(mesh.vertex(iVert), &mesh.targetPointsQueyResults[iVert]);
		//	});
	}
}

//void GAIA::VBDClothSimulationFramework::applyGradientDescent(FloatingType step)
//{
//	if (physicsParams().usePreconditioner)
//	{
//		preconditionGradients();
//
//	}
//	for (IdType iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
//	{
//		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//		mesh.positions() -= mesh.gradientPreconditioned * step;
//	}
//}

void GAIA::VBDClothSimulationFramework::revertGradientDescent(FloatingType step)
{
	for (IdType iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		mesh.positions() += mesh.gradient * step;
	}
}

void GAIA::VBDClothSimulationFramework::printStats(const ConvergenceStats& energyStats, FloatingType energyAll)
{
	std::cout << "Energy at iter " << iIter << ": "
		<< energyAll;

	if (baseTriMeshesForSimulation.size() == 1)
	{
		std::cout << " | inertia energy : " << energyStats.energy_inertia[0]
			<< " | energy_StVK : " << energyStats.energy_elastic_StVK[0]
			<< " | energy_bending : " << energyStats.energy_elastic_bending[0]
			<< " | registration : " << energyStats.energy_registration[0]
			<< " | collision : " << energyStats.energy_collision[0]
			<< " | avg grad norm : " << energyStats.avgGradNorm[0]
			<< " | step size : " << energyStats.stepSize
			<< "\n";
	}
	else
	{
		for (size_t iCloth = 0; iCloth < baseTriMeshesForSimulation.size(); iCloth++)
		{
			std::cout << "Energy for mesh" << iCloth << ": inertia energy: " << energyStats.energy_inertia[iCloth]
				<< " | energy_StVK: " << energyStats.energy_elastic[iCloth]
				<< " | energy_bending : " << energyStats.energy_elastic_bending[iCloth]
				<< " | registration : " << energyStats.energy_registration[iCloth]
				<< " | collision : " << energyStats.energy_collision[iCloth]
				<< " | avg grad norm : " << energyStats.avgGradNorm[iCloth]
				<< " | step size : " << energyStats.stepSize
				<< "\n";
		}
	}
}

void GAIA::VBDClothSimulationFramework::saveStats(const std::string outFile, std::vector<ConvergenceStats> convergenceStats)
{
	nlohmann::json allStats;

	for (size_t iStat = 0; iStat < convergenceStats.size(); iStat++)
	{
		nlohmann::json stat;
		convergenceStats[iStat].toJson(stat);

		allStats.push_back(stat);
	}

	MF::saveJson(outFile, allStats);
}

void GAIA::VBDClothSimulationFramework::updateVelocity()
{
	for (size_t iMesh = 0; iMesh < baseTriMeshesForSimulation.size(); iMesh++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iMesh);
		mesh.velocities = (mesh.positions() - mesh.positionsPrev) / physicsParams().dt;

		cpu_parallel_for(0, mesh.numVertices(), [&](int iVert)
			{
				FloatingType vMag = mesh.velocities.col(iVert).norm();
				FloatingType maxVelocityMagnitude = getObjectParam(iMesh).maxVelocityMagnitude;

				if (getObjectParam(iMesh).hasNoGravZone &&
					maxVelocityMagnitude > 0 &&
					mesh.velocities(GRAVITY_AXIS, iVert) > getObjectParam(iMesh).noGravZoneThreshold)
					// no vel damping in no gravity zone
					// apply the maximum velocity constraint
				{
					if (vMag > maxVelocityMagnitude) {
						mesh.velocities.col(iVert) *= (maxVelocityMagnitude / vMag);
					}
				}
				else if (vMag > 1e-6f) {
					FloatingType vMagNew;
					vMagNew = vMag * getObjectParam(iMesh).exponentialVelDamping - getObjectParam(iMesh).constantVelDamping;
					vMagNew = vMagNew > 1e-6f ? vMagNew : 0.f;
					mesh.velocities.col(iVert) *= vMagNew / vMag;
				}
				else
				{
					mesh.velocities.col(iVert) *= 0.f;
				}
			});

		for (size_t iFixedPoint = 0; iFixedPoint < baseTriMeshesForSimulation[iMesh]->pObjectParams->fixedPoints.size(); iFixedPoint++)
		{
			mesh.velocities.col(mesh.pObjectParams->fixedPoints[iFixedPoint]).setZero();
		}

	}

}


bool GAIA::VBDClothSimulationFramework::determineConvergence(const ConvergenceStats& convergenceStat,
	FloatingType energyAll_prev, FloatingType avgGradNorm_prev, FloatingType currentStepSize)
{
	switch (physicsParams().convergenceType)
	{
	case 0:
		// avgGradNorm + avgGradNormChange
		if (convergenceStat.avgGradNorm_allClothes < physicsParams().convergenceAvgNormThres)
		{
			std::cout << "Convergence Reason: gradient norm smaller than threshold.\n";
			return true;
		}

		else if (abs(avgGradNorm_prev - convergenceStat.avgGradNorm_allClothes) * physicsParams().stepSize / currentStepSize
			< physicsParams().convergenceAvgNormChangeThres)
		{
			std::cout << "Convergence Reason: gradient's change is smaller than threshold.\n";
			return true;
		}
		break;

	case 1:
		// energy
		if (abs(energyAll_prev - convergenceStat.energy_allClothes) < physicsParams().convergenceEnergyChangeThres)
		{
			return true;
		}
		break;
	default:
		break;
	}
	return false;
}

bool GAIA::VBDClothSimulationFramework::writeSimulationParameters(nlohmann::json& outPhysicsParams)
{
	BasePhysicFramework::writeSimulationParameters(outPhysicsParams);
	pClothContactDetectorParameters->toJson(outPhysicsParams["ClothContactDetectorParams"]);
	return false;
}

//void GAIA::VBDClothSimulationFramework::detectIntersections()
//{
//	// bool GAIA::TriTriIntersection::Intersect(IdType fVId1, IdType fVId2, 
//	// const IdType * fVIds1, const IdType * fVIds2, const TVerticesMat & verts1, 
//	// const TVerticesMat & verts2)
//
//	TriTriIntersection tritriIntersection;
//	for (size_t iMesh1 = 0; iMesh1 < baseTriMeshesForSimulation.size(); iMesh1++)
//	{
//		TriMeshFEM::Ptr pMesh1 = baseTriMeshesForSimulation[iMesh1].get();
//		for (size_t iFace1 = 0; iFace1 < pMesh1->numFaces(); iFace1++)
//		{
//			for (size_t iMesh2 = iMesh1; iMesh2 < baseTriMeshesForSimulation.size(); iMesh2++)
//			{
//				TriMeshFEM::Ptr pMesh2 = baseTriMeshesForSimulation[iMesh2].get();
//				size_t startFace = iMesh1 == iMesh2 ? iFace1: 0;
//
//				for (size_t iFace2 = startFace; iFace2 < pMesh2->numFaces(); iFace2++)
//				{
//					if (tritriIntersection.Intersect(iFace1, iFace2, pMesh1->facePos.col(iFace1).data(),
//						pMesh2->facePos.col(iFace2).data(), pMesh1->positions().data(), pMesh2->positions().data())) {
//						std::cout << "Warning: intersection detected between face " << iFace1 << " of mesh " << iMesh1
//							<< " and face " << iFace2 << " of mesh " << iMesh2 << std::endl;
//					}
//				}
//			}
//		}
//	}
//}


void GAIA::VBDClothSimulationFramework::collisionDetection()
{
	if (!physicsParams().handleCollision) {
		return;
	}
	if (iIter == 0)
	{
		pClothContactDetector->updateBVH(RTC_BUILD_QUALITY_LOW);

	}
	else
	{
		pClothContactDetector->updateBVH(RTC_BUILD_QUALITY_REFIT);

	}
	for (int iCloth = 0; iCloth < numClothes(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		cpu_parallel_for(0, mesh.numVertices(), [&](int iV) {
			pClothContactDetector->contactQueryVF(iCloth, iV, &getVFContactResult(iCloth, iV));
			});

		cpu_parallel_for(0, mesh.numEdges(), [&](int iE) {
			pClothContactDetector->contactQueryEE(iCloth, iE, &getEEContactResult(iCloth, iE));
			});
	}

	contactInfoUpdated = true;
}

void GAIA::VBDClothSimulationFramework::applyDeformers()
{
	for (size_t iDeformer = 0; iDeformer < deformers.size(); iDeformer++)
	{
		(*deformers[iDeformer])(*this, curTime, frameId, substep, iIter, physicsParams().dt);
	}
}

size_t GAIA::VBDClothSimulationFramework::numClothes()
{
	return baseTriMeshesForSimulation.size();
}


//FloatingType GAIA::VBDClothSimulationFramework::backTracingLineSearchVBD(FloatingType E0,
//	FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters, FloatingType minStepSizeGD, ConvergenceStats& stats)
//{
//	FloatingType m = 0.f;
//
//	for (int iCloth = 0; iCloth < numClothes(); iCloth++)
//	{
//		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//		orgPosistions[iCloth] = mesh.positions();
//		m += mesh.gradientPreconditioned.squaredNorm();
//	}
//
//	FloatingType bestAlpha = 0.f;
//	FloatingType bestEnergy = E0;
//
//	ConvergenceStats statsTemp(numClothes());
//
//	for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
//	{
//		if (iLineSearchIter != 0)
//		{
//			for (int iCloth = 0; iCloth < numClothes(); iCloth++)
//			{
//				// revert the previous line search step
//				VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//				mesh.positions() = orgPosistions[iCloth];
//			}
//		}
//
//		applyGradientDescent(alpha);
//
//		FloatingType e = evaluateEnergy(statsTemp);
//
//		if (e < bestEnergy)
//		{
//			bestAlpha = alpha;
//			bestEnergy = e;
//			stats = statsTemp;
//		}
//
//		// first Wolfe condition 
//		if (e < E0 - alpha * c * m)
//		{
//			break;
//		}
//		else
//		{
//			alpha = alpha * tau;
//		}
//	}
//
//	if (bestAlpha < minStepSizeGD)
//	{
//		std::cout << "Warning: line search failed, use the smallest step size:" << minStepSizeGD << "\n";
//		bestAlpha = minStepSizeGD;
//
//		for (int iCloth = 0; iCloth < numClothes(); iCloth++)
//		{
//			// revert the previous line search step
//			VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
//			mesh.positions() = orgPosistions[iCloth];
//		}
//		applyGradientDescent(bestAlpha);
//		FloatingType e = evaluateEnergy(stats);
//
//	}
//
//	stats.stepSize = bestAlpha;
//
//	return bestAlpha;
//}

void GAIA::VBDClothSimulationFramework::saveIntermediateResults()
{
	if ((iIter % physicsParams().debugOutputInterval))
	{
		return;
	}

	std::string outputFolderFrame = getDebugFolder() + "/Frame_" + MF::STR::padToNew(std::to_string(frameId), 4, '0') + "_Substep_" + MF::STR::padToNew(std::to_string(substep), 4, '0');
	MF::IO::createFolder(outputFolderFrame);

	std::string outName = outputFolderFrame + +"/iter_" + MF::STR::padToNew(std::to_string(iIter), 5, '0');

	std::string outNameShape = outName + ".ply";
	writeAllToPLY(outNameShape.c_str(), baseTriMeshesForSimulation, basePhysicsParams->saveAllModelsTogether);

	nlohmann::json contactPoints;
	std::vector<std::array<float, 6>> VFContactPointPairs;
	std::vector<std::array<float, 6>> EEContactPointPairs;

	for (int iCloth = 0; iCloth < numClothes(); iCloth++)
	{
		VBDBaseTriMesh& mesh = getSimulatedMesh(iCloth);
		for (size_t iV = 0; iV < mesh.numVertices(); iV++)
		{
			auto& contactResultVF = getVFContactResult(iCloth, iV);
			if (contactResultVF.found) {
				for (size_t iContact = 0; iContact < contactResultVF.contactPts.size(); iContact++)
				{
					VFContactPointPairs.push_back({
						mesh.vertex(iV)(0), mesh.vertex(iV)(1), mesh.vertex(iV)(2),
						contactResultVF.contactPts[iContact].contactPoint(0), contactResultVF.contactPts[iContact].contactPoint(1), contactResultVF.contactPts[iContact].contactPoint(2)
						});
				}
			}
		}

		for (size_t iE = 0; iE < mesh.numEdges(); iE++)
		{
			auto& contactResultEE = getEEContactResult(iCloth, iE);
			if (contactResultEE.found) {
				for (size_t iContact = 0; iContact < contactResultEE.contactPts.size(); iContact++)
				{
					const auto& contactPt = contactResultEE.contactPts[iContact];
					Vec3 C1 = mesh.edgeLerp(iE, contactPt.mu1);
					auto& meshContacting = getMesh(contactPt.contactMeshId2);
					Vec3 C2 = meshContacting.edgeLerp(contactPt.contactEdgeId2, contactPt.mu2);

					EEContactPointPairs.push_back({
						C1(0), C1(1), C1(2),
						C2(0), C2(1), C2(2)
						});
				}
			}
		}

	}
	contactPoints["VFContactPointPairs"] = VFContactPointPairs;
	contactPoints["EEContactPointPairs"] = EEContactPointPairs;
	std::string outNameContacts = outName + ".json";
	MF::saveJson(outNameContacts, contactPoints);

}

void GAIA::VBDClothSimulationFramework::loadTargetMesh(int iFrame)
{
	// load targetMesh
	if (iFrame < physicsParams().targetFiles.size())
	{
		const std::string& targetFile = physicsParams().targetFiles[iFrame];
		targetMeshes.push_back(std::make_shared<TriMeshFEM>());

		// no need to call the initialization function, because it's not used for simulation
		targetMeshes.back()->loadObj(targetFile);
		targetMeshes.back()->computeTopology();

		pTargetMatcher = std::make_shared<MeshClosestPointQuery>(pTargetMatcherParameter);
		pTargetMatcher->initialize(targetMeshes);
	}

}

void GAIA::VBDClothSimulationFramework::newtonEvaluateElasticForceAndHessian()
{
	auto energyPtr = &pNewtonAssembler->elasticEnergy[0];
	auto forcePtr = &pNewtonAssembler->elasticForce[0];
	auto hessianPtr = &pNewtonAssembler->elasticHessian[0];

	for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numFaces(), [&](int iFace) {
			pMesh->computeElasticityForceHessianForFace(iFace, *(energyPtr + iFace), *(forcePtr + iFace), *(hessianPtr + iFace), physicsParams().psd);
			}
		);
		energyPtr += pMesh->numFaces();
		forcePtr += pMesh->numFaces();
		hessianPtr += pMesh->numFaces();
	}
}

void GAIA::VBDClothSimulationFramework::newtonEvaluateCollisionForceAndHessian()
{
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		//for (IdType iV = 0; iV < pMesh->numVertices(); iV++)
		cpu_parallel_for(0, pMesh->numVertices(), [&](IdType iV) {
			// Collision force
			// vertex collision
			ClothVFContactQueryResult& vfContactResult = getVFContactResult(iMesh, iV);

			pNewtonAssembler->vfCollisionInfos[iMesh][iV].resize(vfContactResult.numContactPoints());
			for (size_t iContact = 0; iContact < vfContactResult.numContactPoints(); iContact++)
			{
				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->vfCollisionInfos[iMesh][iV][iContact];
				// vertex is from the v side, so the vertex order is always 3
				collisionInfo.contacting = accumulateVFContactForceAndHessianAllVerts(iContact, &vfContactResult, true, true,
					collisionInfo.lambda, collisionInfo.d2E_dDdD, collisionInfo.normal, collisionInfo.barys);
			}
			});

		cpu_parallel_for(0, pMesh->numEdges(), [&](IdType iE) {
			ClothEEContactQueryResult eeContactResults = getEEContactResult(iMesh, iE);
			pNewtonAssembler->eeCollisionInfos[iMesh][iE].resize(eeContactResults.contactPts.size());
			for (IdType iEEContact = 0; iEEContact < eeContactResults.contactPts.size(); iEEContact++)
			{
				const EEContactPointInfo& eeContact = eeContactResults.contactPts[iEEContact];

				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->eeCollisionInfos[iMesh][iE][iEEContact];
				collisionInfo.contacting = accumulateEEContactForceAndHessianAllVerts(iEEContact, &eeContactResults, true, true,
					collisionInfo.lambda, collisionInfo.d2E_dDdD, collisionInfo.normal, collisionInfo.barys);

			}
			});

	}
}

GAIA::NFloatingType GAIA::VBDClothSimulationFramework::newtonEvaluateCollisionEnergy()
{
	NFloatingType energy = 0;
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		//for (IdType iV = 0; iV < pMesh->numVertices(); iV++)
		cpu_parallel_for(0, pMesh->numVertices(), [&](IdType iV) {
			// Collision force
			// vertex collision
			ClothVFContactQueryResult& vfContactResult = getVFContactResult(iMesh, iV);

			assert(pNewtonAssembler->vfCollisionInfos[iMesh][iV].size() == vfContactResult.numContactPoints());
			for (size_t iContact = 0; iContact < vfContactResult.numContactPoints(); iContact++)
			{
				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->vfCollisionInfos[iMesh][iV][iContact];
				// vertex is from the v side, so the vertex order is always 3
				collisionInfo.contacting = accumulateVFContactEnergy(iContact, &vfContactResult,
					collisionInfo.energy);
			}
			});

		for (int i = 0; i < pMesh->numVertices(); i++) {
			const auto& collisionInfos = pNewtonAssembler->vfCollisionInfos[iMesh][i];
			for (int j = 0; j < collisionInfos.size(); j++) {
				if (collisionInfos[j].contacting) {
					energy += collisionInfos[j].energy;
				}
			}
		}

		cpu_parallel_for(0, pMesh->numEdges(), [&](IdType iE) {
			ClothEEContactQueryResult eeContactResults = getEEContactResult(iMesh, iE);
			assert(pNewtonAssembler->eeCollisionInfos[iMesh][iE].size() == eeContactResults.contactPts.size());
			for (IdType iEEContact = 0; iEEContact < eeContactResults.contactPts.size(); iEEContact++)
			{
				const EEContactPointInfo& eeContact = eeContactResults.contactPts[iEEContact];

				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->eeCollisionInfos[iMesh][iE][iEEContact];
				collisionInfo.contacting = accumulateEEContactEnergy(iEEContact, &eeContactResults,
					collisionInfo.energy);

			}
			});

		for (int i = 0; i < pMesh->numEdges(); i++) {
			const auto& collisionInfos = pNewtonAssembler->eeCollisionInfos[iMesh][i];
			const auto& eeContactResults = getEEContactResult(iMesh, i);
			for (int j = 0; j < collisionInfos.size(); j++) {
				if (collisionInfos[j].contacting) {
					const EEContactPointInfo& eeContact = eeContactResults.contactPts[j];
					if (isSimulationMesh(eeContact.contactMeshId1) && isSimulationMesh(eeContact.contactMeshId2)) {
						energy += collisionInfos[j].energy * 0.5; // same collision is counted twice
					}
					else {
						energy += collisionInfos[j].energy;
					}
				}
			}
		}
	}
	return energy;
}



void GAIA::VBDClothSimulationFramework::fillNewtonSystem()
{
	newtonEvaluateElasticForceAndHessian();
	if (physicsParams().handleCollision) {
		newtonEvaluateCollisionForceAndHessian();
	}

	pNewtonAssembler->newtonForce.setZero();
	fillNewtonForce();

	//pNewtonAssembler->newtonHessian.setZero();
	auto data_ptr = pNewtonAssembler->newtonHessianElasticity.valuePtr();
	if (!pNewtonAssembler->bendingHessianValues.empty()) {
		assert(pNewtonAssembler->newtonHessianElasticity.nonZeros() == pNewtonAssembler->bendingHessianValues.size());
		memcpy(data_ptr, pNewtonAssembler->bendingHessianValues.data(), pNewtonAssembler->bendingHessianValues.size() * sizeof(NFloatingType));
	}
	else {
		memset(data_ptr, 0, pNewtonAssembler->newtonHessianElasticity.nonZeros() * sizeof(NFloatingType));
	}
	if (physicsParams().handleCollision) {
		data_ptr = pNewtonAssembler->newtonHessianCollision.valuePtr();
		memset(data_ptr, 0, pNewtonAssembler->newtonHessianCollision.nonZeros() * sizeof(NFloatingType));
	}

	fillNewtonHessianDiagonal();
	fillNewtonHessianOffDiagonal();
	if (physicsParams().handleCollision) {
		fillNewtonCollisionForceAndHessianV2();
	}
}

void GAIA::VBDClothSimulationFramework::fillNewtonForce()
{
	NVec9* elasticForcePtr = &pNewtonAssembler->elasticForce[0];
	NFloatingType* forcePtr = pNewtonAssembler->newtonForce.data();
	NVecDynamic x;
	x.resize(pNewtonAssembler->numAllVertices * 3);
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {
			if (!pMesh->fixedMask[iV])
			{
				Eigen::Map<NVec3> force(forcePtr + iV * 3);
				// External force and inertia force
				force = (pMesh->vertexExternalForces.col(iV) + pMesh->vertexMass(iV) * (pMesh->inertia.col(iV) - pMesh->vertex(iV)) * (physicsParams().dtSqrReciprocal)).cast<NFloatingType>();
				// Elastic force
				const size_t numNeiFaces = pMesh->numNeiFaces(iV);
				for (size_t iNeiTet = 0; iNeiTet < numNeiFaces; iNeiTet++) {
					const auto faceId = pMesh->getVertexIthNeiFace(iV, iNeiTet);
					const auto vertexFaceOrder = pMesh->getVertexIthNeiFaceOrder(iV, iNeiTet);
					force += elasticForcePtr[faceId].segment<3>(vertexFaceOrder * 3);
				}
			}
			});
		elasticForcePtr += pMesh->numFaces();
		forcePtr += pMesh->numVertices() * 3;

		x.segment(pNewtonAssembler->meshOffsets[iMesh], pMesh->numVertices() * 3) = Eigen::Map<NVecDynamic>(pMesh->positions().data(), pMesh->numVertices() * 3);
	}
	pNewtonAssembler->newtonForce -= pNewtonAssembler->newtonHessianBending * x;
}

void GAIA::VBDClothSimulationFramework::fillNewtonHessianDiagonal()
{
	NMat9* elasticHessianPtr = &pNewtonAssembler->elasticHessian[0];
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {
			NVec3 boundaryForce = NVec3::Zero();
			NMat3 hessian = NMat3::Zero();
			accumlateBoundaryForceAndHessian(pMesh, iMesh, iV, boundaryForce, hessian);

			pNewtonAssembler->newtonForce.segment<3>(iV * 3) += boundaryForce;

			if (!pMesh->fixedMask[iV])
			{
				hessian = NMat3::Identity() * pMesh->vertexMass(iV) * (physicsParams().dtSqrReciprocal);
				// Elastic hessian
				const size_t numNeiTFaces = pMesh->numNeiFaces(iV);
				for (size_t iNeiFace = 0; iNeiFace < numNeiTFaces; iNeiFace++) {
					const auto faceId = pMesh->getVertexIthNeiFace(iV, iNeiFace);
					const auto vertedFaceOrder = pMesh->getVertexIthNeiFaceOrder(iV, iNeiFace);
					hessian += elasticHessianPtr[faceId].block<3, 3>(vertedFaceOrder * 3, vertedFaceOrder * 3);
				}
			}
			else
			{
				hessian = NMat3::Identity() * 1e7f;
			}
			NFloatingType** hessianPtr = &pNewtonAssembler->diagonalHessianBlockPtrs[iMesh][iV * 9];
			for (size_t iRow = 0; iRow < 3; iRow++)
			{
				for (size_t iCol = 0; iCol < 3; iCol++)
				{
					**hessianPtr = hessian(iRow, iCol);
					hessianPtr++;
				}
			}
			});
		elasticHessianPtr += pMesh->numFaces();
	}
}

void GAIA::VBDClothSimulationFramework::fillNewtonHessianOffDiagonal()
{
	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		TriMeshFEM* pMesh = (TriMeshFEM*)baseTriMeshesForSimulation[iMesh].get();
		cpu_parallel_for(0, pMesh->numEdges(), [&](int iE) {

			const EdgeInfo& edge = pMesh->getEdgeInfo(iE);
			int iV = edge.eV1;
			int jV = edge.eV2;
			if (!pMesh->fixedMask[iV] && !pMesh->fixedMask[jV])
			{
				NMat3 hessian = NMat3::Zero();
				// Elastic hessian
				const auto faceId1 = edge.fId1;
				int corner0 = -1, corner1 = -1;
				pMesh->getEdgeVertexOrderInNeighborFace(iE, 0, corner0, corner1);
				hessian += pNewtonAssembler->elasticHessian[faceId1].block<3, 3>(corner0 * 3, corner1 * 3);

				const auto faceId2 = edge.fId2;
				if (faceId2 != -1)
				{
					corner0 = -1, corner1 = -1;
					pMesh->getEdgeVertexOrderInNeighborFace(iE, 1, corner0, corner1);
					hessian += pNewtonAssembler->elasticHessian[faceId2].block<3, 3>(corner0 * 3, corner1 * 3);
				}

				NFloatingType** hessianPtr = &pNewtonAssembler->offDiagonalHessianBlockPtrs[iMesh][iE * 18];
				for (size_t iRow = 0; iRow < 3; iRow++)
				{
					for (size_t iCol = 0; iCol < 3; iCol++)
					{
						**hessianPtr = hessian(iCol, iRow);
						hessianPtr++;
					}
				}
				for (size_t iRow = 0; iRow < 3; iRow++)
				{
					for (size_t iCol = 0; iCol < 3; iCol++)
					{
						**hessianPtr = hessian(iRow, iCol);
						hessianPtr++;
					}
				}
			}
			});
	}
}

void GAIA::VBDClothSimulationFramework::fillNewtonCollisionForceAndHessian()
{
	NVec9* elasticForcePtr = &pNewtonAssembler->elasticForce[0];
	NFloatingType* forcePtr = pNewtonAssembler->newtonForce.data();
	Vec3 normal;
	Vec4 barys;
	FloatingType lambda, d2E_dDdD;

	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		for (IdType iV = 0; iV < pMesh->numVertices(); iV++)
		{
			// Collision force
			// vertex collision
			Mat3 collisionHessian_;
			Vec3 collisionForce = Vec3::Zero();
			ClothVFContactQueryResult& vfContactResult = getVFContactResult(iMesh, iV);

			for (size_t iContact = 0; iContact < vfContactResult.numContactPoints(); iContact++)
			{
				// vertex is from the v side, so the vertex order is always 3
				bool contacting = accumulateVFContactForceAndHessianAllVerts(iContact, &vfContactResult, true, true, lambda, d2E_dDdD, normal, barys);
				if (contacting)
				{
					const VFContactPointInfo& vfContact = vfContactResult.contactPts[iContact];

					CIdType contactVId = vfContact.contactVertexId;
					CIdType contactVertexSideMeshId = vfContact.contactVertexSideMeshId;

					IdType t0 = -1, t1 = -1, t2 = -1;
					if (isSimulationMesh(vfContact.contactFaceSideMeshId))
					{
						const TriMeshFEM* pMeshFSide = baseTriMeshesForSimulation[vfContact.contactFaceSideMeshId].get();
						t0 = pMeshFSide->facePosVId(vfContact.contactFaceId, 0);
						t1 = pMeshFSide->facePosVId(vfContact.contactFaceId, 1);
						t2 = pMeshFSide->facePosVId(vfContact.contactFaceId, 2);
					}

					// must accord with barys' order
					IdType vs[4] = { t0, t1, t2, contactVId };
					IdType meshIds[4] = { vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId, contactVertexSideMeshId };


					for (IdType iRow = 0; iRow < 4; iRow++)
					{
						if (vs[iRow] < 0)
						{
							continue;
						}
						IdType vertPosRow = pNewtonAssembler->meshOffsets[meshIds[iRow]] + vs[iRow] * 3;
						Vec3 collisionForce = lambda * normal * barys[iRow];
						pNewtonAssembler->newtonForce.segment<3>(vertPosRow) += collisionForce;

						for (IdType iCol = 0; iCol < 4; iCol++)
						{
							if (vs[iCol] < 0)
							{
								continue;
							}
							IdType vertPosCol = pNewtonAssembler->meshOffsets[meshIds[iCol]] + vs[iCol] * 3;
							Mat3 hessianBlock = (d2E_dDdD * barys[iRow] * barys[iCol]) * normal * normal.transpose();
							addedToSparse3x3Block(pNewtonAssembler->newtonHessianCollision, vertPosRow, vertPosCol, hessianBlock);
						}
					}
				}
			}

		}

		for (IdType iE = 0; iE < pMesh->numEdges(); iE++)
		{
			ClothEEContactQueryResult eeContactResults = getEEContactResult(iMesh, iE);
			for (IdType iEEContact = 0; iEEContact < eeContactResults.contactPts.size(); iEEContact++)
			{
				const EEContactPointInfo& eeContact = eeContactResults.contactPts[iEEContact];
				TriMeshFEM* pMesh = baseTriMeshesForSimulation[eeContact.contactMeshId1].get();

				IdType e2v1 = -1, e2v2 = -1;
				if (isSimulationMesh(eeContact.contactMeshId2))
				{
					TriMeshFEM* pMeshOtherSide = baseTriMeshesForSimulation[eeContact.contactMeshId2].get();
					e2v1 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV1;
					e2v2 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV2;
				}

				IdType vs[4] = {
					pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV1, pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV2,
					e2v1, e2v2
				};
				IdType meshIds[4] = { eeContact.contactMeshId1, eeContact.contactMeshId1, eeContact.contactMeshId2, eeContact.contactMeshId2 };

				bool contacting = accumulateEEContactForceAndHessianAllVerts(iEEContact, &eeContactResults, true, true, lambda, d2E_dDdD, normal, barys);

				if (contacting)
				{
					for (IdType iRow = 0; iRow < 2; iRow++)
					{
						IdType vertPosRow = pNewtonAssembler->meshOffsets[meshIds[iRow]] + vs[iRow] * 3;

						Vec3 collisionForce = lambda * normal * barys[iRow];
						pNewtonAssembler->newtonForce.segment<3>(vertPosRow) += collisionForce;
						// we only care about the half of the off-diagonal hessian namely hessian(this side, all side), 
						// the other half is handled by the collision results of other side
						for (IdType iCol = 0; iCol < 4; iCol++)
						{
							if (vs[iCol] >= 0 && vs[iRow] >= 0)
							{
								IdType vertPosCol = pNewtonAssembler->meshOffsets[meshIds[iCol]] + vs[iCol] * 3;
								Mat3 hessian = (d2E_dDdD * barys[iRow] * barys[iCol]) * normal * normal.transpose();

								addedToSparse3x3Block(pNewtonAssembler->newtonHessianCollision, vertPosRow, vertPosCol, hessian);
							}
						}
					}
				}
			}
		}
	}
}

void GAIA::VBDClothSimulationFramework::fillNewtonCollisionForceAndHessianV2()
{
	NVec9* elasticForcePtr = &pNewtonAssembler->elasticForce[0];
	NFloatingType* forcePtr = pNewtonAssembler->newtonForce.data();
	Vec3 normal;
	Vec4 barys;
	FloatingType lambda, d2E_dDdD;

	for (size_t iMesh = 0; iMesh < numSimulationMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		for (IdType iV = 0; iV < pMesh->numVertices(); iV++)
		{
			// Collision force
			// vertex collision
			Mat3 collisionHessian_;
			Vec3 collisionForce = Vec3::Zero();
			ClothVFContactQueryResult& vfContactResult = getVFContactResult(iMesh, iV);

			//if (pMesh->fixedMask[iV])
			//{
			//	continue;
			//}

			for (size_t iContact = 0; iContact < vfContactResult.numContactPoints(); iContact++)
			{
				// vertex is from the v side, so the vertex order is always 3
				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->vfCollisionInfos[iMesh][iV][iContact];
				if (collisionInfo.contacting)
				{
					const VFContactPointInfo& vfContact = vfContactResult.contactPts[iContact];

					CIdType contactVId = vfContact.contactVertexId;
					CIdType contactVertexSideMeshId = vfContact.contactVertexSideMeshId;

					IdType t0 = -1, t1 = -1, t2 = -1;
					bool fixed[4] = { true, true,true, pMesh->fixedMask[contactVId] };
					if (isSimulationMesh(vfContact.contactFaceSideMeshId))
					{
						const TriMeshFEM* pMeshFSide = baseTriMeshesForSimulation[vfContact.contactFaceSideMeshId].get();
						t0 = pMeshFSide->facePosVId(vfContact.contactFaceId, 0);
						t1 = pMeshFSide->facePosVId(vfContact.contactFaceId, 1);
						t2 = pMeshFSide->facePosVId(vfContact.contactFaceId, 2);
						fixed[0] = pMesh->fixedMask[t0];
						fixed[1] = pMesh->fixedMask[t1];
						fixed[2] = pMesh->fixedMask[t2];
					}

					// must accord with barys' order
					IdType vs[4] = { t0, t1, t2, contactVId };
					IdType meshIds[4] = { vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId, contactVertexSideMeshId };

					for (IdType iRow = 0; iRow < 4; iRow++)
					{
						if (fixed[iRow])
						{
							continue;
						}
						IdType vertPosRow = pNewtonAssembler->meshOffsets[meshIds[iRow]] + vs[iRow] * 3;
						Vec3 collisionForce = collisionInfo.lambda * collisionInfo.normal * collisionInfo.barys[iRow];
						pNewtonAssembler->newtonForce.segment<3>(vertPosRow) += collisionForce;

						for (IdType iCol = 0; iCol < 4; iCol++)
						{
							if (fixed[iCol])
							{
								continue;
							}
							IdType vertPosCol = pNewtonAssembler->meshOffsets[meshIds[iCol]] + vs[iCol] * 3;
							Mat3 hessianBlock = (collisionInfo.d2E_dDdD * collisionInfo.barys[iRow] * collisionInfo.barys[iCol]) * collisionInfo.normal * collisionInfo.normal.transpose();
							addedToSparse3x3Block(pNewtonAssembler->newtonHessianCollision, vertPosRow, vertPosCol, hessianBlock);
						}
					}
				}
			}

		}

		for (IdType iE = 0; iE < pMesh->numEdges(); iE++)
		{
			ClothEEContactQueryResult eeContactResults = getEEContactResult(iMesh, iE);
			for (IdType iEEContact = 0; iEEContact < eeContactResults.contactPts.size(); iEEContact++)
			{
				const EEContactPointInfo& eeContact = eeContactResults.contactPts[iEEContact];
				TriMeshFEM* pMesh = baseTriMeshesForSimulation[eeContact.contactMeshId1].get();

				IdType e2v1 = -1, e2v2 = -1;
				IdType ev1 = pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV1, ev2 = pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV2;
				bool fixed[4] = { pMesh->fixedMask[ev1],pMesh->fixedMask[ev2], true, true };
				if (isSimulationMesh(eeContact.contactMeshId2))
				{
					TriMeshFEM* pMeshOtherSide = baseTriMeshesForSimulation[eeContact.contactMeshId2].get();
					e2v1 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV1;
					e2v2 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV2;
					fixed[2] = pMeshOtherSide->fixedMask[e2v1];
					fixed[3] = pMeshOtherSide->fixedMask[e2v2];
				}

				IdType vs[4] = {
					ev1, ev2,
					e2v1, e2v2
				};
				IdType meshIds[4] = { eeContact.contactMeshId1, eeContact.contactMeshId1, eeContact.contactMeshId2, eeContact.contactMeshId2 };

				TriMeshCollisionInfoForNewton& collisionInfo = pNewtonAssembler->eeCollisionInfos[iMesh][iE][iEEContact];
				if (collisionInfo.contacting)
				{
					for (IdType iRow = 0; iRow < 2; iRow++)
					{
						if (fixed[iRow])
						{
							continue;
						}
						IdType vertPosRow = pNewtonAssembler->meshOffsets[meshIds[iRow]] + vs[iRow] * 3;

						Vec3 collisionForce = collisionInfo.lambda * collisionInfo.normal * collisionInfo.barys[iRow];
						pNewtonAssembler->newtonForce.segment<3>(vertPosRow) += collisionForce;
						// we only care about the half of the off-diagonal hessian namely hessian(this side, all side), 
						// the other half is handled by the collision results of other side
						for (IdType iCol = 0; iCol < 4; iCol++)
						{
							if (fixed[iCol])
							{
								continue;
							}
							IdType vertPosCol = pNewtonAssembler->meshOffsets[meshIds[iCol]] + vs[iCol] * 3;
							Mat3 hessian = (collisionInfo.d2E_dDdD * collisionInfo.barys[iRow] * collisionInfo.barys[iCol]) * collisionInfo.normal * collisionInfo.normal.transpose();

							addedToSparse3x3Block(pNewtonAssembler->newtonHessianCollision, vertPosRow, vertPosCol, hessian);
						}
					}
				}
			}
		}

	}
}

void GAIA::VBDClothSimulationFramework::newtonEvaluateElasticEnergy()
{
	auto energyPtr = &pNewtonAssembler->elasticEnergy[0];
	for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numFaces(), [&](int iFace) {
			pMesh->computeElasticityEnergyForFace(iFace, *(energyPtr + iFace));
			});
		energyPtr += pMesh->numFaces();
	}
}

void GAIA::VBDClothSimulationFramework::newtonConservativeStepCulling()
{
	for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		cpu_parallel_for(0, pMesh->numVertices(), [&](int vId) {
			if (!pMesh->fixedMask[vId])
			{
				CFloatingType convservativeBounds = pMesh->vertexConvervativeBounds[vId];
				Vec3 accumulatedDisplacement = pMesh->vertex(vId) - pMesh->positionsAtPrevCollisionDetection.col(vId);
				FloatingType accumulatedDisplacementSize = accumulatedDisplacement.norm();
				if (accumulatedDisplacementSize > convservativeBounds)
					// moved out of the conservative bounds, truncation needed, and collision detection is required
				{
					FloatingType accumulatedDisplacementSizeNew = convservativeBounds;
					accumulatedDisplacement = accumulatedDisplacement * (accumulatedDisplacementSizeNew / accumulatedDisplacementSize);
					pMesh->vertex(vId) = pMesh->positionsAtPrevCollisionDetection.col(vId) + accumulatedDisplacement;
					collisionDetectionRequired = true;
				}
			}
			});
	}
}

bool GAIA::VBDClothSimulationFramework::accumulateVFContactForceAndHessianAllVerts(IdType contactId, ClothVFContactQueryResult* pCollisionResult,
	bool updateContactInfo, bool apply_friction, FloatingType& lambda, FloatingType& d2E_dDdD, Vec3& normal, Vec4& barys)
{
	VFContactPointInfo& contactPt = pCollisionResult->contactPts[contactId];

	FloatingType b;
	CFloatingType k = physicsParams().contactStiffness;
	const auto& n = contactPt.contactPointNormal;

	if (updateContactInfo)
	{
		updateVFContactPointInfo(triMeshesAll, contactPt);
	}

	const TriMeshFEM* pMeshVSide = triMeshesAll[contactPt.contactVertexSideMeshId].get();
	const TriMeshFEM* pMeshFSide = triMeshesAll[contactPt.contactFaceSideMeshId].get();

	CVec3 x = pMeshVSide->vertex(contactPt.contactVertexId);
	CVec3 diff = x - contactPt.contactPoint; // from contact point to vertex
	FloatingType dis = diff.dot(n);

	if (dis < physicsParams().contactRadius + physicsParams().thickness)
	{
		barys = { -contactPt.barycentrics[0], -contactPt.barycentrics[1], -contactPt.barycentrics[2], 1.f };

		// simplified point-plane
		// Penalty force
		FloatingType dEdD, d2E_dDdD_;
		computeContactRepulsiveForce(dis, physicsParams().thickness, dEdD, d2E_dDdD_);
		lambda = -dEdD;
		d2E_dDdD = d2E_dDdD_;
		normal = n;
		return true;
	}
	else
	{
		return false;
	}
}

bool GAIA::VBDClothSimulationFramework::accumulateEEContactForceAndHessianAllVerts(IdType contactId, ClothEEContactQueryResult* pCollisionResult,
	bool updateContactInfo, bool apply_friction, FloatingType& lambda, FloatingType& d2E_dDdD, Vec3& normal, Vec4& barys)
{
	EEContactPointInfo& eeContactPt = pCollisionResult->contactPts[contactId];
	if (updateContactInfo)
	{
		updateEEContactPointInfo(triMeshesAll, eeContactPt);
	}

	const Vec3 x1 = eeContactPt.c1;
	const Vec3 x2 = eeContactPt.c2;

	CVec3 diff = x1 - x2; // from contact point 2 to contact point 1
	CFloatingType dis = diff.norm();

	const TriMeshFEM* pMesh1 = triMeshesAll[eeContactPt.contactMeshId1].get();
	const TriMeshFEM* pMesh2 = triMeshesAll[eeContactPt.contactMeshId2].get();
	const EdgeInfo& edgeInfo1 = pMesh1->pTopology->edgeInfos[eeContactPt.contactEdgeId1];
	const EdgeInfo& edgeInfo2 = pMesh2->pTopology->edgeInfos[eeContactPt.contactEdgeId2];

	CVec3 v1 = pMesh1->vertex(edgeInfo1.eV2) - pMesh1->vertex(edgeInfo1.eV1);
	CVec3 v2 = pMesh2->vertex(edgeInfo2.eV2) - pMesh2->vertex(edgeInfo2.eV1);

	CFloatingType parallelEps = v1.cross(v2).squaredNorm() / (v1.norm() * v2.norm());

	if (dis < physicsParams().contactRadius + physicsParams().thickness
		// && dis > 1e-7f 
		&& parallelEps > CMP_EPSILON
		)
	{
		CVec3 n = diff / dis;

		barys = { 1.f - eeContactPt.mu1, eeContactPt.mu1, -1.f + eeContactPt.mu2, -eeContactPt.mu2 };
		CFloatingType k = physicsParams().contactStiffness;
		CFloatingType penetrationDepth = physicsParams().contactRadius - dis;

		FloatingType dEdD, d2E_dDdD_;
		computeContactRepulsiveForce(dis, physicsParams().thickness, dEdD, d2E_dDdD_);
		lambda = -dEdD;
		d2E_dDdD = d2E_dDdD_;
		normal = n;
		return true;
	}
	else
	{
		return false;
	}
}

bool GAIA::VBDClothSimulationFramework::accumulateVFContactEnergy(IdType contactId, ClothVFContactQueryResult* pCollisionResult,
	FloatingType& energy)
{
	VFContactPointInfo& contactPt = pCollisionResult->contactPts[contactId];

	FloatingType b;
	CFloatingType k = physicsParams().contactStiffness;
	const auto& n = contactPt.contactPointNormal;

	const TriMeshFEM* pMeshVSide = triMeshesAll[contactPt.contactVertexSideMeshId].get();
	const TriMeshFEM* pMeshFSide = triMeshesAll[contactPt.contactFaceSideMeshId].get();

	CVec3 x = pMeshVSide->vertex(contactPt.contactVertexId);
	CVec3 diff = x - contactPt.contactPoint; // from contact point to vertex
	FloatingType dis = diff.dot(n);

	if (dis < physicsParams().contactRadius + physicsParams().thickness)
	{
		computeContactEnergy(dis, physicsParams().thickness, energy);
		return true;
	}
	else
	{
		return false;
	}
}

bool GAIA::VBDClothSimulationFramework::accumulateEEContactEnergy(IdType contactId, ClothEEContactQueryResult* pCollisionResult,
	FloatingType& energy)
{
	EEContactPointInfo& eeContactPt = pCollisionResult->contactPts[contactId];
	const Vec3 x1 = eeContactPt.c1;
	const Vec3 x2 = eeContactPt.c2;

	CVec3 diff = x1 - x2; // from contact point 2 to contact point 1
	CFloatingType dis = diff.norm();

	const TriMeshFEM* pMesh1 = triMeshesAll[eeContactPt.contactMeshId1].get();
	const TriMeshFEM* pMesh2 = triMeshesAll[eeContactPt.contactMeshId2].get();
	const EdgeInfo& edgeInfo1 = pMesh1->pTopology->edgeInfos[eeContactPt.contactEdgeId1];
	const EdgeInfo& edgeInfo2 = pMesh2->pTopology->edgeInfos[eeContactPt.contactEdgeId2];

	CVec3 v1 = pMesh1->vertex(edgeInfo1.eV2) - pMesh1->vertex(edgeInfo1.eV1);
	CVec3 v2 = pMesh2->vertex(edgeInfo2.eV2) - pMesh2->vertex(edgeInfo2.eV1);

	CFloatingType parallelEps = v1.cross(v2).squaredNorm() / (v1.norm() * v2.norm());

	if (dis < physicsParams().contactRadius + physicsParams().thickness
		// && dis > 1e-7f 
		&& parallelEps > CMP_EPSILON
		)
	{
		CFloatingType k = physicsParams().contactStiffness;
		CFloatingType penetrationDepth = physicsParams().contactRadius - dis;

		computeContactEnergy(dis, physicsParams().thickness, energy);
		return true;
	}
	else
	{
		return false;
	}
}

GAIA::NFloatingType GAIA::VBDClothSimulationFramework::newtonEvaluateMeritEnergy(NFloatingType& eInertia, NFloatingType& eElastic, NFloatingType& eBending, NFloatingType& eContact, bool elasticReady)
{
	if (!elasticReady)newtonEvaluateElasticEnergy();
	eElastic = 0.f;
	for (auto& e : pNewtonAssembler->elasticEnergy)
	{
		eElastic += e;
	}
	eInertia = 0;
	eBending = 0;
	eContact = 0;
	NVecDynamic x;
	x.resize(pNewtonAssembler->numAllVertices * 3);
	for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	{
		VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
		eInertia += pMesh->evaluatedMeritEnergyInertia();
		x.segment(pNewtonAssembler->meshOffsets[iMesh], pMesh->numVertices() * 3) = Eigen::Map<NVecDynamic>(pMesh->positions().data(), pMesh->numVertices() * 3);
	}
	eBending = 0.5 * x.transpose() * pNewtonAssembler->newtonHessianBending * x;

	if (physicsParams().handleCollision) {
		eContact = newtonEvaluateCollisionEnergy();
	}
	return eElastic + eInertia + eBending + eContact;
}

GAIA::NFloatingType GAIA::VBDClothSimulationFramework::newtonLineSearch(const VecDynamic& dx, NFloatingType E0, FloatingType alpha,
	FloatingType c, FloatingType tau, int maxNumIters, FloatingType& stepSizeOut)
{
	FloatingType m = dx.squaredNorm();

	std::vector<TVerticesMat> orgPos{};
	orgPos.reserve(numTriMeshes());
	for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	{
		orgPos.push_back(getSimulatedMesh(iMesh).positions());
	}

	NFloatingType e = E0;
	NFloatingType eInertia = 0;
	NFloatingType eElastic = 0;
	NFloatingType eBending = 0;
	NFloatingType eContact = 0;
	bool satisfied = false;
	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
		std::cout << "Initial Energy: " << E0 << std::endl;
		});
	for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
	{
		int offset = 0;
		for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
		{
			VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
			pMesh->positions() = orgPos[iMesh] + alpha * Eigen::Map<const TVerticesMat>(dx.data() + offset, 3, pMesh->numVertices());
			offset += pMesh->numVertices() * 3;
		}
		e = newtonEvaluateMeritEnergy(eInertia, eElastic, eBending, eContact);
		debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
			std::cout << "alpha: " << alpha << ", energy: " << e << ", inertia: " << eInertia << ", elastic: " << eElastic << ", bending: " << eBending << ", contact: " << eContact << std::endl;
			});

		// first Wolfe condition 
		if (e < E0 - alpha * c * m)
		{
			// std::cout << "step size for vertex " << vId << ": " << alpha << "\n";
			satisfied = true;
			break;
		}
		else
		{
			alpha = alpha * tau;
		}
	}
	stepSizeOut = alpha;
	return e;
	//int offset = 0;
	//for (size_t iMesh = 0; iMesh < numTriMeshes(); iMesh++)
	//{
	//	VBDBaseTriMesh* pMesh = &getSimulatedMesh(iMesh);
	//	pMesh->positions() = orgPos[iMesh];
	//	offset += pMesh->numVertices() * 3;
	//}
	//stepSizeOut = 0;
	//return E0;
}
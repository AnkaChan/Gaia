
#include "../CollisionDetector/DiscreteCollisionDetector.h"
#include "../CollisionDetector/ContinuousCollisionDetector.h"

#include "VBDPhysics.h"
#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>
#include <MeshFrame/Utility/Time.h>


#include "../Timer/Timer.h"
#include "../Timer/RunningTimeStatistics.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"

#include "VBDPhysicsCompute.h"
#include "VBD_GeneralCompute.h"

#include "../../3rdParty/CuMatrix/CuMatrix/MatrixOps/CuMatrixVis.h"
#include "../Utility/checksum.h"

#include "VBDDeformer.h"

#include "../BVH/MortonCode.h"

// #define APPLY_LOCAL_LINE_SEARCH

using namespace GAIA;

//GAIA::ObjectParamsVBD& GAIA::ObjectParamsEBDList::getObjectParam(int iObj)
//{
//	return *objectParams[iObj];
//}
//
//bool GAIA::ObjectParamsEBDList::fromJson(nlohmann::json& objectParam)
//{
//	nlohmann::json modelsInfo = objectParam["Models"];
//
//	for (nlohmann::json& model : modelsInfo) {
//		std::string materialName;
//		parseJsonParameters(model, "materialName", materialName);
//
//		if (materialName == "NeoHookean")
//		{
//			objectParams.push_back(std::make_shared<ObjectParametersEBDNeoHookean>());
//			objectParams.back()->fromJson(model);
//			objectParams.back()->materialName = materialName;
//			// tMeshes.push_back(std::make_shared<PBDTetMeshNeoHookean>());
//		}
//		else if (materialName == "MassSpring") {
//			// wait to be fullfilled
//			
//		}
//		else
//		{
//			std::cout << "Warning!!! Material name: " << materialName << " not recognized! Skipping this object!\n";
//			/*std::cout << "Using default material NeoHookean instead!\n";
//			objectParams.push_back(std::make_shared<ObjectParamsPBDNeoHookean>());*/
//		}
//	}
//
//	return true;
//}
//
//bool GAIA::ObjectParamsEBDList::toJson(nlohmann::json& objectParam)
//{
//	for (size_t iObj = 0; iObj < objectParams.size(); iObj++)
//	{
//		nlohmann::json obj;
//		objectParams[iObj]->toJson(obj);
//		objectParam["Models"].push_back(obj);
//	}
//	return true;
//}
//
//bool GAIA::EBDPhysicsAllParameters::fromJson(nlohmann::json& physicsJsonParams)
//{
//	pPhysicsParams.fromJson(physicsJsonParams["PhysicsParams"]);
//	collisionParams.fromJson(physicsJsonParams["CollisionParams"]);
//	return true;
//}
//
//bool GAIA::EBDPhysicsAllParameters::toJson(nlohmann::json& physicsJsonParams)
//{
//	pPhysicsParams.toJson(physicsJsonParams["PhysicsParams"]);
//	collisionParams.toJson(physicsJsonParams["CollisionParams"]);
//	return true;
//}


void GAIA::VBDPhysics::loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder)
{
	BasePhysicFramework::loadRunningparameters(inModelInputFile, inParameterFile, outFolder);

	physicsParamsVBD = std::static_pointer_cast<VBDPhysicsParameters>(basePhysicsParams);


}

inline IdType findSmallestParallelGroup(std::vector<IdType>& availableColors, const std::vector<std::vector<IdType>>& vertexParallelGroups,
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

struct TetInfo
{
	IdType meshId;
	IdType tetId;
	IdType vertedOrderInTet;
	IdType vertexId;
	bool used = false;
};

size_t getNumSharedVerts(const TetInfo & tet1, const TetInfo& tet2, const std::vector<VBDBaseTetMesh::SharedPtr>& tMeshes) 
{
	if (tet1.meshId != tet2.meshId)
	{
		return 0;
	}
	else
	{
		VBDBaseTetMesh::Ptr pTM = tMeshes[tet1.meshId].get();

		IdType* tetVIds1 = pTM->tetVIds().col(tet1.tetId).data();
		IdType* tetVIds2 = pTM->tetVIds().col(tet2.tetId).data();

		size_t numSharedVerts = 0;

		for (size_t i = 0; i < 4; i++)
		{
			for (size_t j = 0; j < 4; j++)
			{
				if (tetVIds1[i] == tetVIds2[i]) {
					++numSharedVerts;
					break;
				}
			}
		}

		return numSharedVerts;
	}
}

IdType getMostAfftinitiveTet(const std::vector<TetInfo>& currentTetColorGroupSorted, const std::vector<TetInfo>& currentTetColorGroup, size_t numTetsToConsider, const std::vector<VBDBaseTetMesh::SharedPtr>& tMeshes) {
	if (currentTetColorGroupSorted.size()==0)
	{
		return 0;
	}

	IdType maxNumOfSharedVerts = -1;
	IdType mostAfftinitiveTet = -1;

	for (IdType iTet = 0; iTet < currentTetColorGroup.size(); iTet++)
	{
		IdType numOfSharedVerts = 0;
		const TetInfo& tet1 = currentTetColorGroup[iTet];

		if (tet1.used)
		{
			continue;
		}

		numTetsToConsider = std::min(numTetsToConsider, currentTetColorGroupSorted.size());

		for (IdType i = numTetsToConsider; i >= 1 ; i--)

		{
			const TetInfo& tet2 = currentTetColorGroupSorted[currentTetColorGroupSorted.size()-i];

			numOfSharedVerts += getNumSharedVerts(tet1, tet2, tMeshes);
			
		}

		if (numOfSharedVerts > maxNumOfSharedVerts)
		{
			maxNumOfSharedVerts = numOfSharedVerts;
			mostAfftinitiveTet = iTet;
		}
	}

	// std::cout << "Most Afftinitive Tet: " << currentTetColorGroup[mostAfftinitiveTet].tetId << " number of shared verts: " << maxNumOfSharedVerts << std::endl;

	assert(mostAfftinitiveTet != -1);

	return mostAfftinitiveTet;
}


void GAIA::VBDPhysics::initialize()
{
	BasePhysicFramework::initialize();

	for (size_t iMesh = 0; iMesh < basetetMeshes.size(); iMesh++)
	{
		tMeshes.push_back(getTetMeshSharedPtrAs<VBDBaseTetMesh>(iMesh));
		objectParamsVBD.push_back(std::static_pointer_cast<ObjectParamsVBD>(objectParamsList->objectParams[iMesh]));
		debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
			MaterialForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			InertiaForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			boundaryFrictionForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			boundaryCollisionForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			frictionForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			collisionForce.push_back(TVerticesMat::Zero(3, tMeshes[iMesh]->numVertices()));
			});
	}

	// initialize collision results
	collisionResultsAll.resize(numTetMeshes());
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		collisionResultsAll[iMesh].resize(basetetMeshes[iMesh]->surfaceVIds().size());
	}

	// sortVertexColorGroupsByMortonCode();

	// generate parallelization groups
	size_t numberOfParallelGroups = 0;
	for (size_t iMesh = 0; iMesh < basetetMeshes.size(); iMesh++)
	{
		if (numberOfParallelGroups < basetetMeshes[iMesh]->verticesColoringCategories().size())
		{
			numberOfParallelGroups = basetetMeshes[iMesh]->verticesColoringCategories().size();
		}
	}

	vertexParallelGroups.resize(numberOfParallelGroups);

	numAllVertices = 0;
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		size_t numColors = pMesh->verticesColoringCategories().size();
		std::vector<IdType> availableColors(numColors, -1);
		numAllVertices += pMesh->numVertices();

		for (int iColor = 0; iColor < numColors; iColor++)
		{
			availableColors[iColor] = iColor;
		}

		for (int iColor = 0; iColor < numColors; iColor++)
		{
			const std::vector<IdType>& currentColorGroup = pMesh->verticesColoringCategories()[iColor];

			int smallestGroupId = findSmallestParallelGroup(availableColors, vertexParallelGroups, true);
			// add the current color group to this parallel group
			for (int iVertex = 0; iVertex < currentColorGroup.size(); iVertex++)
			{
				vertexParallelGroups[smallestGroupId].push_back(iMesh);
				vertexParallelGroups[smallestGroupId].push_back(currentColorGroup[iVertex]);
			}
			// color group from the same mesh cannot be added to this group again
			// std::cout << "Color " << iColor << " of mesh " << iMesh << " was added to parallel group " << smallestGroupId << "\n";
		}
	}

	// assemble tet parallel group
	// preallocate the parallel group collision list
	tetParallelGroups.resize(vertexParallelGroups.size());
	//for (int iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
	std::mutex printMutex;
	cpu_parallel_for(0, vertexParallelGroups.size(), [&](int iGroup)
		{
			const auto& currentVertexColorGroup = vertexParallelGroups[iGroup];
			auto& currentTetColorGroup = tetParallelGroups[iGroup];
			currentTetColorGroup.reserve(currentVertexColorGroup.size());
			{
				std::lock_guard<std::mutex> printLockGuard(printMutex);
				std::cout << "size of vertex parallel group " << iGroup << ": " << currentVertexColorGroup.size() / 2 << "\n";
			}
			
			std::vector<std::vector<TetInfo>> currentTetColorGroupPerMesh;
			currentTetColorGroupPerMesh.resize(tMeshes.size());
			
			// std::vector<TetInfo> currentTetColorGroupTemp;

			for (int iVertex = 0; iVertex < currentVertexColorGroup.size() / 2; iVertex++)
			{
				int meshId = currentVertexColorGroup[iVertex * 2];
				VBDBaseTetMesh::SharedPtr pMesh = tMeshes[meshId];
				int vertexId = currentVertexColorGroup[iVertex * 2 + 1];
				// set the parallel group id of each vertex
				pMesh->vertexParallelGroups(vertexId) = iGroup;

				const size_t numNeiTest = pMesh->getNumVertexNeighborTets(vertexId);

				//std::vector<float4> neiTetCentroids;
				for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
				{
					IdType tetId = pMesh->getVertexNeighborTet(vertexId, iNeiTet);
					int vertedOrderInTet = pMesh->getVertexNeighborTetVertexOrder(vertexId, iNeiTet);

					TetInfo tetInfo;
					tetInfo.meshId = meshId;
					tetInfo.tetId = tetId;
					tetInfo.vertedOrderInTet = vertedOrderInTet;
					tetInfo.vertexId = vertexId;
					
					currentTetColorGroupPerMesh[meshId].push_back(tetInfo);
					//currentTetColorGroupTemp.push_back(tetInfo);

					//IdType* tetVIds = pMesh->tet(tetId).data();
					//Vec3 centroid = (pMesh->vertex(tetVIds[0]) + pMesh->vertex(tetVIds[1])
					//	+ pMesh->vertex(tetVIds[2]) + pMesh->vertex(tetVIds[3])) * 0.25;

					//float4 p_;
					//p_.x = centroid.x();
					//p_.y = centroid.y();
					//p_.z = centroid.z();

					//neiTetCentroids.push_back(p_);
				}

				//// sort each vertex's local neighbor tets based on morton code
				//std::vector<IdType> sortedIndices;
				//mortonSortingCPU(neiTetCentroids, sortedIndices);

				//std::vector<IdType> orderedNeiTets;
				//std::vector<IdType> orderedVertedOrders;
				//for (size_t i = 0; i < numNeiTest; i++)
				//{
				//	orderedNeiTets.push_back(pMesh->getVertexNeighborTet(vertexId, sortedIndices[i]));
				//	orderedVertedOrders.push_back(pMesh->getVertexNeighborTetVertexOrder(vertexId, sortedIndices[i]));
				//}
				//for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
				//{
				//	pMesh->pTopology->vertexNeighborTets(pMesh->pTopology->vertexNeighborTets_infos(2 * vertexId) + iNeiTet) = orderedNeiTets[iNeiTet];
				//	pMesh->pTopology->vertexNeighborTets_vertexOrder(pMesh->pTopology->vertexNeighborTets_infos(2 * vertexId) + iNeiTet) = orderedVertedOrders[iNeiTet];
				//}

			}

			for (int meshId = 0; meshId < currentTetColorGroupPerMesh.size(); meshId++)
			{
				std::vector<float4> ps;
				VBDBaseTetMesh::SharedPtr pMesh = tMeshes[meshId];
				std::vector<TetInfo>& currentTetColorGroupThisMesh = currentTetColorGroupPerMesh[meshId];
				for (int iTet = 0; iTet < currentTetColorGroupThisMesh.size(); ++iTet)
				{
					IdType* tetVIds = pMesh->tet(currentTetColorGroupThisMesh[iTet].tetId).data();

					Vec3 centroid = (pMesh->vertex(tetVIds[0]) + pMesh->vertex(tetVIds[1]) 
						+ pMesh->vertex(tetVIds[2]) + pMesh->vertex(tetVIds[3])) * 0.25;

					float4 p_;
					p_.x = centroid.x();
					p_.y = centroid.y();
					p_.z = centroid.z();

					ps.push_back(p_);
				}

				std::vector<IdType> sortedIndices;
				lbvh::mortonSortingCPU(ps, sortedIndices);

				for (size_t iTet = 0; iTet < sortedIndices.size(); iTet++)
				{
					const TetInfo& tet = currentTetColorGroupThisMesh[sortedIndices[iTet]];
					currentTetColorGroup.push_back(tet.meshId);
					currentTetColorGroup.push_back(tet.tetId);
					currentTetColorGroup.push_back(tet.vertedOrderInTet);
					currentTetColorGroup.push_back(tet.vertexId);
				}
			}

			//std::vector<TetInfo> currentTetColorGroupSorted;

			//currentTetColorGroupSorted.reserve(currentTetColorGroupTemp.size());
			//while (currentTetColorGroupTemp.size()!= currentTetColorGroupSorted.size()) {
			//	IdType mostAffinitiveTetId = getMostAfftinitiveTet(currentTetColorGroupSorted, currentTetColorGroupTemp, NUM_THREADS_TET_SWEEP/2, tMeshes);

			//	currentTetColorGroupSorted.push_back(currentTetColorGroupTemp[mostAffinitiveTetId]);

			//	currentTetColorGroupTemp[mostAffinitiveTetId].used = true;
			//}

			//for (size_t iTet = 0; iTet < currentTetColorGroupSorted.size(); iTet++)
			//{
			//	const TetInfo& tet = currentTetColorGroupSorted[iTet];
			//	currentTetColorGroup.push_back(tet.meshId);
			//	currentTetColorGroup.push_back(tet.tetId);
			//	currentTetColorGroup.push_back(tet.vertedOrderInTet);
			//	currentTetColorGroup.push_back(tet.vertexId);
			//}
				
			{
				std::lock_guard<std::mutex> printLockGuard(printMutex);
				std::cout << "size of tet parallel group " << iGroup << ": " << currentTetColorGroup.size() / 4 << "\n";
			}
		});


	activeColllisionList.initialize(tMeshes, vertexParallelGroups, physicsParams().activeCollisionListPreAllocationRatio);
	if (physicsParams().useGPU)
	{
		initializeGPU();
	}

	if (physicsParams().debugGPU)
	{
		initializeGPU_cpuDebugData();
	}

	// load deformers
	auto deformerParams = physicsJsonParams["Deformers"];

	for (auto deformerParam : deformerParams) {
		deformers.push_back(loadDeformers(*this, deformerParam));
	}

	if (physicsParams().useNewton)
	{
		initializeNewton();
	}

	if (physicsParams().useGDSolver)
	{
		initializeGD();
		if (pLineSearchUtilities == nullptr)
		{
			pLineSearchUtilities = std::make_shared<LineSearchUtilities>();
			pLineSearchUtilities->initialize(*this);
		}
	}

	if (physicsParams().useLineSearch || physicsParams().evaluateConvergence)
	{
		pNewtonAssembler->elasticEnergy.resize(numAllTets);
		if (pLineSearchUtilities == nullptr)
		{
			pLineSearchUtilities = std::make_shared<LineSearchUtilities>();
			pLineSearchUtilities->initialize(*this);
		}
	}
}

void GAIA::VBDPhysics::initializeNewton()
{
	pNewtonAssembler = std::make_shared<TetMeshNewtonAssembler>();

	pNewtonAssembler->diagonalHessianBlockPtrs.resize(tMeshes.size());
	pNewtonAssembler->offDiagonalHessianBlockPtrs.resize(tMeshes.size());
	std::vector<NTriplet> newtonHessianTriplets;
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];

		pNewtonAssembler->diagonalHessianBlockPtrs[iMesh].reserve(pMesh->numVertices() * 9);
		pNewtonAssembler->offDiagonalHessianBlockPtrs[iMesh].reserve(pMesh->numEdges() * 18);

	}
	newtonHessianTriplets.reserve(numAllVertices * 9 + numAllEdges * 18);
	int offset = 0;
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTriplets.emplace_back(vertPos + j, vertPos + i, 1.0);
				}
			}
		}
		for (int iE = 0; iE < pMesh->numEdges(); iE++)
		{
			int vertPosi = offset + pMesh->edges()(0, iE) * 3;
			int vertPosj = offset + pMesh->edges()(1, iE) * 3;
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					newtonHessianTriplets.emplace_back(vertPosi + i, vertPosj + j, 1.0);
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTriplets.emplace_back(vertPosj + j, vertPosi + i, 1.0);
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}
	pNewtonAssembler->newtonHessianAll.resize(numAllVertices * 3, numAllVertices * 3);
	pNewtonAssembler->newtonHessianAll.setFromTriplets(newtonHessianTriplets.begin(), newtonHessianTriplets.end());
	pNewtonAssembler->newtonHessianAll.makeCompressed();

	if (physicsParams().NewtonUseCG)
	{
		pNewtonAssembler->solverCG.compute(pNewtonAssembler->newtonHessianAll);
		pNewtonAssembler->solverCG.setMaxIterations(300);
		pNewtonAssembler->solverCG.setTolerance(1e-7f);
	}
	else {
		pNewtonAssembler->solverDirect.analyzePattern(pNewtonAssembler->newtonHessianAll);
	}

	offset = 0;
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					pNewtonAssembler->diagonalHessianBlockPtrs[iMesh].push_back(&pNewtonAssembler->newtonHessianAll.coeffRef(vertPos + j, vertPos + i));
				}
			}
		}
		for (int iE = 0; iE < pMesh->numEdges(); iE++)
		{
			int vertPosi = offset + pMesh->edges()(0, iE) * 3;
			int vertPosj = offset + pMesh->edges()(1, iE) * 3;
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					pNewtonAssembler->offDiagonalHessianBlockPtrs[iMesh].push_back(&pNewtonAssembler->newtonHessianAll.coeffRef(vertPosi + i, vertPosj + j));
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					pNewtonAssembler->offDiagonalHessianBlockPtrs[iMesh].push_back(&pNewtonAssembler->newtonHessianAll.coeffRef(vertPosj + j, vertPosi + i));
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}
	pNewtonAssembler->newtonForce.resize(numAllVertices * 3);
	pNewtonAssembler->elasticHessian.resize(numAllTets);
	pNewtonAssembler->elasticForce.resize(numAllTets);
}

void GAIA::VBDPhysics::initializeGPU()
{
	cudaStreamCreate(&cudaStream);
	tetMehesGPUBuffer = std::make_shared<ManagedBuffer<VBDBaseTetMeshGPU*>>(tMeshes.size(), true);
	tMeshesGPU.resize(tMeshes.size(), nullptr);

	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		//cpu_parallel_for(0, tMeshes.size(), [&](int iMesh)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		pMesh->initializeGPUMesh();
		tetMehesGPUBuffer->getCPUBuffer()[iMesh] = pMesh->getGPUMesh();
		tMeshesGPU[iMesh] = pMesh->getGPUMesh();
	}
	//);
	tetMehesGPUBuffer->toGPU(false, cudaStream);
	//parallelGroupsHeadGPUBuffer = std::make_shared<ManagedBuffer<int32_t*>>(vertexParallelGroups.size(), true);

	for (size_t i = 0; i < vertexParallelGroups.size(); i++)
	{
		vertexParallelGroupsBuffer.push_back(
			std::make_shared<ManagedBuffer<int32_t>>(vertexParallelGroups[i].size(),
				true,
				vertexParallelGroups[i].data())
		);

		activeCollisionsEachParallelGroupBuffer.push_back(
			std::make_shared<ManagedBuffer<int32_t>>(activeColllisionList.activeCollisionsEachParallelGroup[i].size(),
				true,
				activeColllisionList.activeCollisionsEachParallelGroup[i].data())
		);

		//parallelGroupsHeadGPUBuffer->getCPUBuffer()[i] = vertexParallelGroupsBuffer.back()->getGPUBuffer();
		vertexParallelGroupsBuffer.back()->toGPU(false, cudaStream);
		vertexParallelGroupHeadsGPU.push_back(vertexParallelGroupsBuffer.back()->getGPUBuffer());

		tetParallelGroupsBuffer.push_back(
			std::make_shared<ManagedBuffer<int32_t>>(tetParallelGroups[i].size(),
				true,
				tetParallelGroups[i].data())
		);
		tetParallelGroupsBuffer.back()->toGPU(false, cudaStream);
		tetParallelGroupHeadsGPU.push_back(tetParallelGroupsBuffer.back()->getGPUBuffer());
	}

	vertexAllParallelGroupsBuffer = std::make_shared<ManagedBuffer<int32_t>>(2 * numAllVertices, true);

	size_t vertexAccumulate = 0;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		for (size_t iV = 0; iV < pMesh->numVertices(); iV++)
		{
			vertexAllParallelGroupsBuffer->getCPUBuffer()[vertexAccumulate * 2] = iMesh;
			vertexAllParallelGroupsBuffer->getCPUBuffer()[vertexAccumulate * 2 + 1] = iV;
			++vertexAccumulate;
		}
	}
	assert(vertexAccumulate == numAllVertices);
	vertexAllParallelGroupsBuffer->toGPU();

	vbdPhysicsDataGPUCPUBuffer.boundaryCollisionStiffness = physicsParams().boundaryCollisionStiffness;
	vbdPhysicsDataGPUCPUBuffer.boundaryFrictionDynamic = physicsParams().boundaryFrictionDynamic;
	vbdPhysicsDataGPUCPUBuffer.boundaryFrictionEpsV = physicsParams().boundaryFrictionEpsV;

	vbdPhysicsDataGPUCPUBuffer.collisionStiffness = physicsParams().collisionStiffness;
	vbdPhysicsDataGPUCPUBuffer.collisionAirDistance = physicsParams().collisionAirDistance;

	//vbdPhysicsDataGPUCPUBuffer.collisionOffHeight = physicsParams().collisionOffHeight;

	vbdPhysicsDataGPUCPUBuffer.dt = physicsParams().dt;
	vbdPhysicsDataGPUCPUBuffer.dtSqrReciprocal = physicsParams().dtSqrReciprocal;
	vbdPhysicsDataGPUCPUBuffer.stepSize = physicsParams().stepSize;
	vbdPhysicsDataGPUCPUBuffer.stepSizeGD = physicsParams().stepSizeGD;
	vbdPhysicsDataGPUCPUBuffer.worldBounds[0] = physicsParams().worldBounds(0, 0);
	vbdPhysicsDataGPUCPUBuffer.worldBounds[1] = physicsParams().worldBounds(1, 0);
	vbdPhysicsDataGPUCPUBuffer.worldBounds[2] = physicsParams().worldBounds(2, 0);
	vbdPhysicsDataGPUCPUBuffer.worldBounds[3] = physicsParams().worldBounds(0, 1);
	vbdPhysicsDataGPUCPUBuffer.worldBounds[4] = physicsParams().worldBounds(1, 1);
	vbdPhysicsDataGPUCPUBuffer.worldBounds[5] = physicsParams().worldBounds(2, 1);
	vbdPhysicsDataGPUCPUBuffer.tetMeshes = tetMehesGPUBuffer->getGPUBuffer();
	vbdPhysicsDataGPUCPUBuffer.numMeshes = numTetMeshes();
	vbdPhysicsDataGPUCPUBuffer.useAccelerator = physicsParams().useAccelerator;
	vbdPhysicsDataGPUCPUBuffer.useBlockJacobi = physicsParams().GDSolverUseBlockJacobi;

	pPhysicsDataGPUBuffer = std::make_shared<DeviceClassBuffer<VBDPhysicsDataGPU>>();
	pPhysicsDataGPUBuffer->fromCPU(&vbdPhysicsDataGPUCPUBuffer);
	CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}

void GAIA::VBDPhysics::initializeGPU_cpuDebugData()
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		tetMehesGPUBuffer_forCPUDebug.push_back(pMesh->getGPUMesh_forCPUDebug());
	}
	vbdPhysicsDataGPU_forCPUDebug.tetMeshes = tetMehesGPUBuffer_forCPUDebug.data();
	vbdPhysicsDataGPU_forCPUDebug.numMeshes = numTetMeshes();

	vbdPhysicsDataGPU_forCPUDebug.boundaryCollisionStiffness = physicsParams().boundaryCollisionStiffness;
	vbdPhysicsDataGPU_forCPUDebug.boundaryFrictionDynamic = physicsParams().boundaryFrictionDynamic;
	vbdPhysicsDataGPU_forCPUDebug.boundaryFrictionEpsV = physicsParams().boundaryFrictionEpsV;

	//vbdPhysicsDataGPU_forCPUDebug.collisionOffHeight = physicsParams().collisionOffHeight;

	vbdPhysicsDataGPU_forCPUDebug.collisionStiffness = physicsParams().collisionStiffness;
	vbdPhysicsDataGPU_forCPUDebug.collisionAirDistance = physicsParams().collisionAirDistance;
	vbdPhysicsDataGPU_forCPUDebug.dt = physicsParams().dt;
	vbdPhysicsDataGPU_forCPUDebug.dtSqrReciprocal = physicsParams().dtSqrReciprocal;
	vbdPhysicsDataGPU_forCPUDebug.stepSize = physicsParams().stepSize;
	vbdPhysicsDataGPU_forCPUDebug.stepSizeGD = physicsParams().stepSizeGD;
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[0] = physicsParams().worldBounds(0, 0);
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[1] = physicsParams().worldBounds(1, 0);
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[2] = physicsParams().worldBounds(2, 0);
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[3] = physicsParams().worldBounds(0, 1);
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[4] = physicsParams().worldBounds(1, 1);
	vbdPhysicsDataGPU_forCPUDebug.worldBounds[5] = physicsParams().worldBounds(2, 1);
	vbdPhysicsDataGPU_forCPUDebug.useAccelerator = physicsParams().useAccelerator;
	vbdPhysicsDataGPU_forCPUDebug.useBlockJacobi = physicsParams().GDSolverUseBlockJacobi;
}

void GAIA::VBDPhysics::enableModels()
{

	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pTetMesh = tMeshes[iMesh];
		if (pTetMesh->pObjectParams->frameToAppear <= frameId) {
			bool changed = false;
			if (false == pTetMesh->activeForCollision)
			{
				pTetMesh->activeForCollision = true;
				changed = true;
			}
			if (false == pTetMesh->activeForMaterialSolve)
			{
				pTetMesh->activeForMaterialSolve = true;
				changed = true;
			}
			if (changed && physicsParams().useGPU)
			{
				pTetMesh->setGPUMeshActiveness(true, true);
			}
		}
	}
}

void GAIA::VBDPhysics::initializeGD()
{
	pGDSolverUtilities = std::make_shared<GD_SolverUtilities>();
	pGDSolverUtilities->initialize(*this);
}

void GAIA::VBDPhysics::sortVertexColorGroupsByMortonCode()
{
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh::SharedPtr pMesh = tMeshes[iMesh];
		size_t numColors = pMesh->verticesColoringCategories().size();

		for (int iColor = 0; iColor < numColors; iColor++)
		{
			std::vector<IdType>& currentColorGroup = pMesh->verticesColoringCategories()[iColor];

			std::vector<float4> ps;
			for (IdType vId : currentColorGroup)
			{
				Vec3 p = pMesh->vertex(vId);
				float4 p_;
				p_.x = p.x();
				p_.y = p.y();
				p_.z = p.z();

				ps.push_back(p_);
			}

			std::vector<IdType> sortedIndices;

			lbvh::mortonSortingCPU(ps, sortedIndices);

			std::vector<IdType> currentColorGroupNew(currentColorGroup.size());

			for (size_t i = 0; i < sortedIndices.size(); i++)
			{
				currentColorGroupNew[i] = currentColorGroup[sortedIndices[i]];
			}

			currentColorGroup = std::move(currentColorGroupNew);
		}
	}
}



TetMeshFEM::SharedPtr GAIA::VBDPhysics::initializeMaterial(ObjectParams::SharedPtr objParam, std::shared_ptr<TetMeshMF> pTMeshMF,
	BasePhysicsParams::SharedPtr physicsParaemters)
{
#ifdef KEEP_MESHFRAME_MESHES
	tMeshesMF[iMesh] = pTM_MF;
#endif // KEEP_MESHFRAME_MESHES
	TetMeshFEM::SharedPtr pBaseMesh;
	switch (objParam->materialType)
	{
	case NeoHookean:
	{
		VBDTetMeshNeoHookean::SharedPtr pTetMeshNeoHookean = std::make_shared<VBDTetMeshNeoHookean>();
		pBaseMesh = pTetMeshNeoHookean;
		pTetMeshNeoHookean->initialize(objParam, pTMeshMF, physicsParamsVBD, this);
		//GPUTMeshes[iMesh] = pTetMeshNeoHookean->getTetMeshGPU();
		break;
	}
	case MassSpring:
	{
		//EBDTetMesh_MassSpring::SharedPtr pTetMeshMassSpring = std::make_shared<EBDTetMesh_MassSpring>();
		//pBaseMesh = pTetMeshMassSpring;
		//pTetMeshMassSpring->initialize(objParam, pTMeshMF, physicsParamsVBD, this);

		//GPUTMeshes[iMesh] = pTetMeshMassSpring->getTetMeshGPU();

		std::cout << "Warning!!! Material name: " << objParam->materialName << " not implemented yet! Skipping this object!\n";
		assert(false);
	}
	break;
	default:
		break;
	}

	return pBaseMesh;
}

void GAIA::VBDPhysics::timeStepEqualize()
{
	BasePhysicFramework::timeStepEqualize();

	for (size_t iMesh = 0; iMesh < objectParamsList->size(); iMesh++)
	{
		// make parameters step invariant
		objectParamsList->getObjectParamAs<ObjectParamsVBD>(iMesh).exponentialVelDamping =
			std::pow(objectParamsList->getObjectParamAs<ObjectParamsVBD>(iMesh).exponentialVelDamping, physicsParams().dt);
		objectParamsList->getObjectParamAs<ObjectParamsVBD>(iMesh).constantVelDamping =
			objectParamsList->getObjectParamAs<ObjectParamsVBD>(iMesh).constantVelDamping * physicsParams().dt;

	}
}

void GAIA::VBDPhysics::runStep()
{
	if (physicsParams().useNewton)
	{
		runStepNewton();
	}
	else if (physicsParams().useGDSolver)
	{
		runStepGPU_GD();
	}
	else if (physicsParams().useGPU)
	{
		runStepGPU();
		//runStepGPU_allInOneSweep();
	}
	else {
		switch (physicsParams().collisionSolutionType)
		{
		case 0:
			runStep_serialCollisionHandling();
			break;
		case 1:
			runStep_hybridCollisionHandling();
			break;
		default:
			break;
		}
	}
}

void GAIA::VBDPhysics::runStepGPU()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});
		dcd();

		TICK(timeCsmpInitialStep);
		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
			if (tMeshes[iMesh]->activeForMaterialSolve)
			{
				tMeshes[iMesh]->evaluateExternalForce();
				tMeshes[iMesh]->applyInitialStep();
			}
			});
		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
		ccd();

		prepareCollisionDataGPU();

		applyDeformers();

		syncAllToGPU(false);

		if (physicsParams().evaluateConvergence)
		{
			iIter = -1;
			evaluateConvergenceGPU();
		}

		FloatingType acceleratorOmega = 1.f;
		if (physicsParams().useAccelerator)
		{
			recordInitialPositionForAccelerator(false);
		}

		TICK(timeCsmpMaterialSolve);
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			if (physicsParams().intermediateCollisionIterations > 0
				&& !(iIter % physicsParams().intermediateCollisionIterations)
				&& iIter)
			{
				intermediateCollisionDetection();
			}
			if (physicsParams().useAccelerator)
			{
				acceleratorOmega = getAcceleratorOmega(iIter + 1, physicsParams().acceleratorPho, acceleratorOmega);
			}

			if (!graphCreated)
			{
				recordVBDSolveGraph();
			}

			if (physicsParams().useAccelerator)
			{
				recordPrevIterPositionsAccelerator(false);
			}
			cudaGraphLaunch(VBDSolveInstance, cudaStream);

			//for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
			//{
			//	// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
			//	// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);

			//	VBDSolveParallelGroup_tetSweep(getVBDPhysicsDataGPU(), tetParallelGroupHeadsGPU[iGroup],
			//		sizeTetParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);
			//	//CHECK_CUDA_ERROR(cudaDeviceSynchronize());

			//	size_t numActiveCollisions = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
			//	if (numActiveCollisions)
			//	{
			//		VBDSolveParallelGroup_collisionSweep(getVBDPhysicsDataGPU(), numActiveCollisions,
			//			activeCollisionsEachParallelGroupBuffer[iGroup]->getGPUBuffer(), physicsParams().numThreadsVBDSolve, cudaStream);
			//	}
			//	//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
			//	VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
			//		sizeVertexParallelGroup(iGroup), acceleratorOmega, physicsParams().numThreadsVBDSolve, cudaStream);

			//	//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
			//}

			//if (iIter && (iIter % 10 == 0))
			//{
			//	CHECK_CUDA_ERROR(cudaDeviceSynchronize());

			//}

			//if (physicsParams().useAccelerator && acceleratorOmega != 1.f)
			//{
			//	// copy accelerated positions to vertPos
			//	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
			//	{
			//		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
			//		CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshShared->vertPosBuffer->getGPUBuffer(), pMesh->pTetMeshSharedBase->positionsNewBuffer->getGPUBuffer(),
			//			pMesh->pTetMeshShared->positionsNewBuffer->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
			//	}
			//}

			if (physicsParams().evaluateConvergence && (iIter % physicsParams().evaluationSteps == 0))
			{
				evaluateConvergenceGPU();
			}

		} // iteration

		// copy the CPU data before syncing
		//TVerticesMat vertFromCPU = tMeshes[0]->positions();

		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;

		TICK(timeCsmpUpdateVelocity);
		//updateVelocities();
		updateVelocitiesGPU();
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);

		syncAllToCPU(true);
		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);

	} // substep
}

//void GAIA::VBDPhysics::runStepGPU_allInOneSweep()
//{
//	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
//	{
//		curTime += physicsParams().dt;
//		debugOperation(DEBUG_LVL_DEBUG, [&]() {
//			std::cout << "Substep step: " << substep << std::endl;
//			});
//		dcd();
//
//		TICK(timeCsmpInitialStep);
//		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
//			if (tMeshes[iMesh]->activeForMaterialSolve)
//			{
//				tMeshes[iMesh]->evaluateExternalForce();
//				tMeshes[iMesh]->applyInitialStep();
//			}
//			});
//		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
//		ccd();
//
//		prepareCollisionDataGPU();
//
//		applyDeformers();
//
//		syncAllToGPU(false);
//
//		if (physicsParams().evaluateConvergence)
//		{
//			iIter = -1;
//			evaluateConvergenceGPU();
//		}
//
//		FloatingType acceleratorOmega = 1.f;
//		if (physicsParams().useAccelerator)
//		{
//			recordInitialPositionForAccelerator(false);
//		}
//
//		TICK(timeCsmpMaterialSolve);
//		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
//		{
//			if (physicsParams().intermediateCollisionIterations > 0
//				&& !(iIter % physicsParams().intermediateCollisionIterations)
//				&& iIter)
//			{
//				intermediateCollisionDetection();
//			}
//			if (physicsParams().useAccelerator)
//			{
//				acceleratorOmega = getAcceleratorOmega(iIter + 1, physicsParams().acceleratorPho, acceleratorOmega);
//			}
//
//			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
//			{
//				// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
//				// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);
//
//				VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
//					sizeVertexParallelGroup(iGroup), acceleratorOmega, physicsParams().numThreadsVBDSolve, cudaStream);
//
//				//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//			}
//
//			if (physicsParams().useAccelerator && acceleratorOmega != 1.f)
//			{
//				// copy accelerated positions to vertPos
//				for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
//				{
//					VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
//					CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshShared->vertPosBuffer->getGPUBuffer(), pMesh->pTetMeshSharedBase->positionsNewBuffer->getGPUBuffer(),
//						pMesh->pTetMeshShared->positionsNewBuffer->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
//				}
//			}
//
//			if (physicsParams().evaluateConvergence && iIter % physicsParams().evaluationSteps == 0)
//			{
//				evaluateConvergenceGPU();
//			}
//
//		} // iteration
//
//		// copy the CPU data before syncing
//		//TVerticesMat vertFromCPU = tMeshes[0]->positions();
//
//		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;
//
//		TICK(timeCsmpUpdateVelocity);
//		//updateVelocities();
//		updateVelocitiesGPU();
//		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);
//
//		syncAllToCPU(true);
//		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);
//
//	} // substep
//}

//void GAIA::VBDPhysics::runStepGPU_acceleratedGS()
//{
//	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
//	{
//		curTime += physicsParams().dt;
//		debugOperation(DEBUG_LVL_DEBUG, [&]() {
//			std::cout << "Substep step: " << substep << std::endl;
//			});
//		dcd();
//
//		TICK(timeCsmpInitialStep);
//		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
//			if (tMeshes[iMesh]->activeForMaterialSolve)
//			{
//				tMeshes[iMesh]->evaluateExternalForce();
//				tMeshes[iMesh]->applyInitialStep();
//			}
//			});
//		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
//		ccd();
//
//		prepareCollisionDataGPU();
//
//		applyDeformers();
//
//		syncAllToGPU(false);
//
//		if (physicsParams().evaluateConvergence)
//		{
//			iIter = -1;
//			evaluateConvergenceGPU();
//		}
//
//		FloatingType acceleratorOmega = 1.f;
//		if (physicsParams().useAccelerator)
//		{
//			recordInitialPositionForAccelerator(false);
//		}
//
//		TICK(timeCsmpMaterialSolve);
//		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
//		{
//			if (physicsParams().intermediateCollisionIterations > 0
//				&& !(iIter % physicsParams().intermediateCollisionIterations)
//				&& iIter)
//			{
//				intermediateCollisionDetection();
//			}
//			if (physicsParams().useAccelerator)
//			{
//				acceleratorOmega = getAcceleratorOmega(iIter + 1, physicsParams().acceleratorPho, acceleratorOmega);
//			}
//
//			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
//			{
//				// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
//				// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);
//
//				VBDSolveParallelGroup_tetSweep(getVBDPhysicsDataGPU(), tetParallelGroupHeadsGPU[iGroup],
//					sizeTetParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);
//				//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//
//				size_t numActiveCollisions = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
//				if (numActiveCollisions)
//				{
//					VBDSolveParallelGroup_collisionSweep(getVBDPhysicsDataGPU(), numActiveCollisions,
//						activeCollisionsEachParallelGroupBuffer[iGroup]->getGPUBuffer(), physicsParams().numThreadsVBDSolve, cudaStream);
//				}
//				//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//				if (physicsParams().useAccelerator)
//				{
//					VBDSolveParallelGroup_vertexSweepAcceleratedGS(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
//						sizeVertexParallelGroup(iGroup), acceleratorOmega, physicsParams().numThreadsVBDSolve, cudaStream);
//				}
//				else {
//					VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
//						sizeVertexParallelGroup(iGroup), cudaStream);
//				}
//
//
//				//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//			}
//
//			if (physicsParams().evaluateConvergence && iIter % physicsParams().evaluationSteps == 0)
//			{
//				evaluateConvergenceGPU();
//			}
//
//		} // iteration
//
//		// copy the CPU data before syncing
//		//TVerticesMat vertFromCPU = tMeshes[0]->positions();
//
//		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;
//
//		TICK(timeCsmpUpdateVelocity);
//		//updateVelocities();
//		updateVelocitiesGPU();
//		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);
//
//		syncAllToCPU(true);
//		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);
//
//	} // substep
//}

void GAIA::VBDPhysics::runStepGPU_GD()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});
		//dcd();

		TICK(timeCsmpInitialStep);
		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
			if (tMeshes[iMesh]->activeForMaterialSolve)
			{
				tMeshes[iMesh]->evaluateExternalForce();
				tMeshes[iMesh]->applyInitialStep();
			}
			});
		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
		//ccd();

		//prepareCollisionDataGPU();

		applyDeformers();

		syncAllToGPU(false);

		if (physicsParams().evaluateConvergence)
		{
			iIter = -1;
			evaluateConvergenceGPU();
		}

		TICK(timeCsmpMaterialSolve);

		// compute initial energy
		pLineSearchUtilities->ePrev = evaluateMeritEnergyGPU(pLineSearchUtilities->eInertiaPrev, pLineSearchUtilities->eElasticPrev);

		if (pGDSolverUtilities->stepSizePrevStep < physicsParams().stepSizeGD)
		{
			pGDSolverUtilities->stepSizePrevStep /= physicsParams().backtracingLineSearchTau;
		}

		FloatingType acceleratorOmega = 1.f;
		if (physicsParams().useAccelerator)
		{
			recordInitialPositionForAccelerator(false);
		}
		GDBackupPositions(false, acceleratorOmega);

		int lineSearchStartIter = -1;

		// syncAllToCPUVertPosOnly(true);
		// saveDebugState("Initialization_", true);

		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			if (physicsParams().useAccelerator)
			{
				acceleratorOmega = getAcceleratorOmega(iIter + 1, physicsParams().acceleratorPho, acceleratorOmega);
			}
			//syncAllToCPUVertPosOnly(true);
			//saveDebugState("beforeIter_retry_" + std::to_string(retry) + "_", true);

			if (physicsParams().GDSolverUseBlockJacobi)
			{
				GDSolveParallelGroup_BlockJacobi_allInOneSweepGPU(getVBDPhysicsDataGPU(), vertexAllParallelGroupsBuffer->getGPUBuffer(), numAllVertices,
					physicsParams().numThreadsVBDSolve, cudaStream);
			}
			else {

				GDSolveParallelGroup_allInOneSweepGPU(getVBDPhysicsDataGPU(), vertexAllParallelGroupsBuffer->getGPUBuffer(), numAllVertices,
					physicsParams().numThreadsVBDSolve, cudaStream);
			}
			GDSolveParallelGroup_updatePositionSweepGPU(getVBDPhysicsDataGPU(), pGDSolverUtilities->stepSizePrevStep, acceleratorOmega, vertexAllParallelGroupsBuffer->getGPUBuffer(), numAllVertices,
				physicsParams().numThreadsVBDSolve, cudaStream);
			//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
			//syncAllToCPUVertPosOnly(true);
			//saveDebugState("retry_" + std::to_string(retry) + "_", true);

			if (physicsParams().useLineSearch)
			{
				if (physicsParams().lineSearchGapIter == 1 || (iIter && iIter % physicsParams().lineSearchGapIter == 0)) {

					pLineSearchUtilities->e = evaluateMeritEnergyGPU(pLineSearchUtilities->eInertia, pLineSearchUtilities->eElastic);

					if (std::isnan(pLineSearchUtilities->e) || pLineSearchUtilities->e > pLineSearchUtilities->ePrev)
					{
						pGDSolverUtilities->stepSizePrevStep *= physicsParams().backtracingLineSearchTau;
						GDRevertToBackupPositions(false, acceleratorOmega);
						iIter = lineSearchStartIter;

					}
					else
					{
						GDBackupPositions(false, acceleratorOmega);
						pLineSearchUtilities->recordPrevEnergy();
						lineSearchStartIter = iIter;
					}

				}
				// //std::cout << "Line search at iIter: " << iIter
				// //	<< " | e: " << e << " | eInertia: " << eInertia << " : eElastic: " << eElastic << std::endl;
				// syncAllToCPUVertPosOnly(true);
				// // compute line search
				// cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
				// 	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
				// 	pGDSolverUtilities->dxs[iMesh] = pMesh->vertices() - pGDSolverUtilities->previousPositionsForLineSearch[iMesh];
				// });

				// GDLineSearch(pGDSolverUtilities, ePrev, eInertiaPrev, eElasticPrev, physicsParams().backtracingLineSearchAlpha, physicsParams().backtracingLineSearchC,
				// 	physicsParams().backtracingLineSearchTau, physicsParams().backtracingLineSearchMaxIters);
				// syncAllToGPUVertPosOnly(false);
			}

			if (physicsParams().evaluateConvergence && iIter % physicsParams().evaluationSteps == 0)
				// if this got called before 
			{
				evaluateConvergenceGPU();
			}

		} // iteration

		// copy the CPU data before syncing
		//TVerticesMat vertFromCPU = tMeshes[0]->positions();

		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;

		TICK(timeCsmpUpdateVelocity);
		//updateVelocities();
		updateVelocitiesGPU();
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);

		syncAllToCPU(true);
		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);

	} // substep
}

//void GAIA::VBDPhysics::runStepGPU_debugOnCPU()
//{
//
//	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
//	{
//		debugOperation(DEBUG_LVL_DEBUG, [&]() {
//			std::cout << "Substep step: " << substep << std::endl;
//			});
//		dcd();
//
//		TICK(timeCsmpInitialStep);
//		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
//			if (tMeshes[iMesh]->activeForMaterialSolve)
//			{
//				tMeshes[iMesh]->evaluateExternalForce();
//				tMeshes[iMesh]->applyInitialStep();
//			}
//			});
//		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
//		ccd();
//
//		prepareCollisionDataCPUAndGPU();
//		testGPUCollisionHandlingCode();
//
//		TICK(timeCsmpMaterialSolve);
//		syncAllToGPU(false);
//		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
//		{
//			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
//			{
//				// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
//				// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);
//
//				VBDSolveParallelGroup_tetSweep(getVBDPhysicsDataGPU(), tetParallelGroupHeadsGPU[iGroup],
//					sizeTetParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);
//
//				size_t numActiveCollisions = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
//				if (numActiveCollisions)
//				{
//					VBDSolveParallelGroup_collisionSweep(getVBDPhysicsDataGPU(), numActiveCollisions,
//						activeCollisionsEachParallelGroupBuffer[iGroup]->getGPUBuffer(), physicsParams().numThreadsVBDSolve, cudaStream);
//				}
//				//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//
//				VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
//					sizeVertexParallelGroup(iGroup), cudaStream);
//
//				// CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//			}
//
//		} // iteration
//
//		// copy the CPU data before syncing
//		//TVerticesMat vertFromCPU = tMeshes[0]->positions();
//
//		syncAllToCPUVertPosOnly(true);
//		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);
//
//		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;
//
//		TICK(timeCsmpUpdateVelocity);
//		updateVelocities();
//		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);
//	} // substep
//}

//void GAIA::VBDPhysics::runStepGPUNoCollision()
//{
//	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
//	{
//		debugOperation(DEBUG_LVL_DEBUG, [&]() {
//			std::cout << "Substep step: " << substep << std::endl;
//			});
//
//		TICK(timeCsmpInitialStep);
//		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
//			if (tMeshes[iMesh]->activeForMaterialSolve)
//			{
//				tMeshes[iMesh]->evaluateExternalForce();
//				tMeshes[iMesh]->applyInitialStep();
//			}
//			});
//		TOCK_STRUCT(timeStatistics(), timeCsmpInitialStep);
//
//		TICK(timeCsmpMaterialSolve);
//
//		syncAllToGPU(false);
//
//		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
//		{
//			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
//			{
//				// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
//				// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);
//
//				VBDSolveParallelGroup_tetSweep(getVBDPhysicsDataGPU(), tetParallelGroupHeadsGPU[iGroup],
//					sizeTetParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);
//
//				VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
//					sizeVertexParallelGroup(iGroup), cudaStream);
//
//				// CHECK_CUDA_ERROR(cudaDeviceSynchronize());
//			}
//
//		} // iteration
//
//		// copy the CPU data before syncing
//		//TVerticesMat vertFromCPU = tMeshes[0]->positions();
//
//		syncAllToCPUVertPosOnly(true);
//		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);
//
//		//std::cout << vertFromCPU - tMeshes[0]->positions() << std::endl;
//
//		TICK(timeCsmpUpdateVelocity);
//		updateVelocities();
//		TOCK_STRUCT(timeStatistics(), timeCsmpUpdateVelocity);
//	} // substep
//}

void GAIA::VBDPhysics::runStepNewton()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});

		// dcd();

		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
			if (tMeshes[iMesh]->activeForMaterialSolve)
			{
				tMeshes[iMesh]->evaluateExternalForce();
				tMeshes[iMesh]->applyInitialStep();
			}
			});
		applyDeformers();

		// ccd();
		// prepareCollisionDataCPU();

		if (physicsParams().evaluateConvergence)
		{
			iIter = -1;
			evaluateConvergence();
		}

		FloatingType stepSize = physicsParams().backtracingLineSearchAlpha;
		TICK(timeCsmpMaterialSolve);

		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			bool apply_friction = iIter >= physicsParams().frictionStartIter;
			apply_friction = true;
			debugOperation(DEBUG_LVL_DEBUG, [&]() {
				std::cout << "iIter: " << iIter << std::endl;
				});
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE_2, std::bind(&VBDPhysics::outputPosVel, this));
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE_2, std::bind(&VBDPhysics::clearForces, this));
			
			computeElasticForceHessian();
			fillNewtonSystem();
			NVecDynamic Ndx ;

			if (physicsParams().NewtonUseCG)
			{
				pNewtonAssembler->solverCG.compute(pNewtonAssembler->newtonHessianAll);
				Ndx = pNewtonAssembler->solverCG.solve(pNewtonAssembler->newtonForce);

				if (pNewtonAssembler->solverCG.info() != Eigen::Success)
				{
					std::cerr << "CG solve failed. Error: " << pNewtonAssembler->solverCG.info() << std::endl;
					//std::exit(-1);
				}
			}
			else
			{
				pNewtonAssembler->solverDirect.factorize(pNewtonAssembler->newtonHessianAll);
				Ndx = pNewtonAssembler->solverDirect.solve(pNewtonAssembler->newtonForce);

				if (pNewtonAssembler->solverDirect.info() != Eigen::Success)
				{
					std::cerr << "Factorization failed. Error: " << pNewtonAssembler->solverDirect.info() << std::endl;
					std::exit(-1);
				}
			}

			VecDynamic dx = Ndx.cast<FloatingType>();

			debugOperation(DEBUG_LVL_DEBUG, [&]() {
				NCFloatingType averageForceNorm = Eigen::Map<NTVerticesMat>(pNewtonAssembler->newtonForce.data(), 3, numAllVertices).colwise().norm().mean();
				std::cout << "averageForceNorm: " << averageForceNorm << std::endl;
				});
			if (dx.hasNaN()) {
				std::cerr << "Newton system has NaNs" << std::endl;
				std::exit(-1);
			}
			if (physicsParams().useLineSearch) {
				if (iIter == 0) {
					NFloatingType eInertia{};
					NFloatingType eElastic{};
					// Elastic is already computed during force and hessian evaluation, set elasticReady to true
					pNewtonAssembler->newtonEnergy = evaluateMeritEnergy(eInertia, eElastic, true);
				}
				FloatingType stepSizeNew;
				pNewtonAssembler->newtonEnergy = newtonLineSearch(dx, pNewtonAssembler->newtonEnergy, stepSize, physicsParams().backtracingLineSearchC,
					physicsParams().backtracingLineSearchTau, physicsParams().backtracingLineSearchMaxIters, stepSizeNew);
				stepSize = physicsParams().backtracingLineSearchAlpha;
			}
			else {
				updatePositions(dx);
			}
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE_2, std::bind(&VBDPhysics::outputForces, this));
			//if (physicsParams().intermediateCollisionIterations > 0 && iIter % physicsParams().intermediateCollisionIterations == physicsParams().intermediateCollisionIterations - 1) {
			//	intermediateCollisionDetection();

			if (physicsParams().evaluateConvergence && iIter % physicsParams().evaluationSteps == 0)
			{
				evaluateConvergence();
			}

		} // iteration
		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);



		updateVelocities();
	} // substep
}

void GAIA::VBDPhysics::recordVBDSolveGraph()
{
	cudaStreamBeginCapture(cudaStream, cudaStreamCaptureModeGlobal);
	for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
	{
		// VBDSolveParallelGroup_allInOneSweep(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup], 
		// 	sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve);

		//VBDSolveParallelGroup_tetSweep(getVBDPhysicsDataGPU(), tetParallelGroupHeadsGPU[iGroup],
		//	sizeTetParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);
		//CHECK_CUDA_ERROR(cudaDeviceSynchronize());

		VBDSolveParallelGroup_vertexSweep_2hierarchiesGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
			sizeVertexParallelGroup(iGroup), cudaStream);

		VBDSolveParallelGroup_updateVertexPositionGPU(getVBDPhysicsDataGPU(), vertexParallelGroupHeadsGPU[iGroup],
			sizeVertexParallelGroup(iGroup), physicsParams().numThreadsVBDSolve, cudaStream);

		//CHECK_CUDA_ERROR(cudaDeviceSynchronize());
	}
	cudaStreamEndCapture(cudaStream, &VBDSolveGraph);
	cudaGraphInstantiate(&VBDSolveInstance, VBDSolveGraph, NULL, NULL, 0);
	graphCreated = true;

}

void GAIA::VBDPhysics::applyDeformers()
{
	for (size_t iDeformer = 0; iDeformer < deformers.size(); iDeformer++)
	{
		(*deformers[iDeformer])(*this, curTime, frameId, substep, iIter, physicsParams().dt);
	}
}

void GAIA::VBDPhysics::runStep_serialCollisionHandling()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});

		dcd();

		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
			if (tMeshes[iMesh]->activeForMaterialSolve)
			{
				tMeshes[iMesh]->evaluateExternalForce();
				tMeshes[iMesh]->applyInitialStep();
			}
			});
		applyDeformers();

		ccd();

		prepareCollisionDataCPU();

		TICK(timeCsmpMaterialSolve);
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			bool apply_friction = iIter >= physicsParams().frictionStartIter;
			solveCollisionsSequentially();
			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
			{
				const std::vector<IdType>& parallelGroup = vertexParallelGroups[iGroup];

				size_t numVertices = parallelGroup.size() / 2;
				cpu_parallel_for(0, numVertices, [&](int iV) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];

					VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
					if (!pMesh->fixedMask[vId] && !pMesh->activeCollisionMask[vId] && pMesh->activeForMaterialSolve)
						//if (!pMesh->fixedMask[vId])
					{
						//pMesh->VBDStep(vId);
						VBDStepWithCollision(pMesh, iMesh, vId, apply_friction);
					}
					});
			}
		} // iteration
		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);

		updateVelocities();
	} // substep
}

void GAIA::VBDPhysics::runStep_hybridCollisionHandling()
{
	for (substep = 0; substep < physicsParams().numSubsteps; substep++)
	{
		curTime += physicsParams().dt;
		debugOperation(DEBUG_LVL_DEBUG, [&]() {
			std::cout << "Substep step: " << substep << std::endl;
			});

		dcd();

		cpu_parallel_for(0, numTetMeshes(), [&](int iMesh) {
			if (tMeshes[iMesh]->activeForMaterialSolve)
			{
				tMeshes[iMesh]->evaluateExternalForce();
				tMeshes[iMesh]->applyInitialStep();
			}
			});
		applyDeformers();

		ccd();
		prepareCollisionDataCPU();


		TICK(timeCsmpMaterialSolve);
		for (iIter = 0; iIter < physicsParams().iterations; iIter++)
		{
			bool apply_friction = iIter >= physicsParams().frictionStartIter;
			// apply_friction = true;
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
				std::cout << "iIter: " << iIter << std::endl;
				});
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, std::bind(&VBDPhysics::outputPosVel, this));
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, std::bind(&VBDPhysics::clearForces, this));
			for (size_t iGroup = 0; iGroup < vertexParallelGroups.size(); iGroup++)
			{
				const std::vector<IdType>& parallelGroup = vertexParallelGroups[iGroup];

				size_t numVertices = parallelGroup.size() / 2;
				size_t numCollisionParallelGroup = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
				// update collision info for parallel group
				cpu_parallel_for(0, numCollisionParallelGroup, [&](int iCollision) {
					IdType iMesh = activeColllisionList.activeCollisionsEachParallelGroup[iGroup][iCollision * 2];
					int vId = activeColllisionList.activeCollisionsEachParallelGroup[iGroup][2 * iCollision + 1];
					updateCollisionInfo(getCollisionDetectionResultFromTetMeshId(iMesh, vId));

					});

				cpu_parallel_for(0, numVertices, [&](int iV) {
					IdType iMesh = parallelGroup[iV * 2];
					int vId = parallelGroup[2 * iV + 1];

					VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
					if (!pMesh->fixedMask[vId] && pMesh->activeForMaterialSolve)
						//if (!pMesh->fixedMask[vId])
					{
						//pMesh->VBDStep(vId);
						VBDStepWithCollision(pMesh, iMesh, vId, apply_friction);
					}
					});
				// VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[0].get();
				//int sum = checksum(reinterpret_cast<int*>(pMesh->mVertPos.data()), sizeof(FloatingType) / sizeof(int) * pMesh->mVertPos.size());
				//std::cout << "substep: " << substep << "iteration: " << iteration << "iGroup: " << iGroup << ", checksum: " << sum << std::endl;
			}
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, std::bind(&VBDPhysics::outputForces, this));
			if (physicsParams().intermediateCollisionIterations > 0 && iIter % physicsParams().intermediateCollisionIterations == physicsParams().intermediateCollisionIterations - 1) {
				intermediateCollisionDetection();
			}
		} // iteration
		TOCK_STRUCT(timeStatistics(), timeCsmpMaterialSolve);

		updateVelocities();
	} // substep
}

void GAIA::VBDPhysics::prepareCollisionDataCPU()
{
	// record the activate collisions
	activeColllisionList.clear();

	int numActiveDCD = 0;
	int numActiveCCD = 0;

	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++) {
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		pTetMesh->activeCollisionMask.setZero();
		std::vector<VBDCollisionDetectionResult>& dcdResults = collisionResultsAll[iMesh];

		// todo: make this parallel using 
		for (int iSurfaceV = 0; iSurfaceV < pTetMesh->surfaceVIds().size(); ++iSurfaceV)
		{
			VBDCollisionDetectionResult& collisionResult = dcdResults[iSurfaceV];
			int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
			if (collisionResult.numIntersections())
			{
				activeColllisionList.addToActiveCollisionList(collisionResult);
				//std::cout << "collision mask of " << surfaceVId << " set to true\n";

				if (collisionResult.fromCCD)
				{
					numActiveCCD++;
				}
				else {
					numActiveDCD++;
				}

				for (size_t iIntersection = 0; iIntersection < collisionResult.numIntersections(); iIntersection++)
				{
					const CollidingPointInfo& collidingPt = collisionResult.collidingPts[iIntersection];
					if (collidingPt.shortestPathFound)
					{
						pTetMesh->activeCollisionMask(surfaceVIdTetMesh) = true;

						int meshId_intersecting = collidingPt.intersectedMeshId;
						int surfaceFaceId = collidingPt.closestSurfaceFaceId;
						VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();

						//std::cout << "surfaceFaceId " << surfaceFaceId << "\n";

						for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
						{
							IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);
							//std::cout << "collision mask of " << faceVId << " set to true\n";
							CollisionRelationList& collisionRelation = getCollisionRelationList(meshId_intersecting, faceVId);
							pTetMesh_intersecting->activeCollisionMask(faceVId) = true;
							collisionRelation.emplace_back();
							collisionRelation.back().meshId = iMesh;
							collisionRelation.back().surfaceVertexId = iSurfaceV;
							collisionRelation.back().collisionId = iIntersection;
							collisionRelation.back().collisionType = 0;
							collisionRelation.back().collisionVertexOrder = iSurfaceFaceV;
						}
					}
				}
			}
		}
	}
	// std::cout << "-------------------------------------\n";
	//for (size_t iCol = 0; iCol < activeColllisionList.activeCollisions.size(); iCol++)
	//{
	//	int meshId = activeColllisionList.activeCollisions[iCol].first;
	//	int vertexTetMeshId = activeColllisionList.activeCollisions[iCol].second;
	//	std::cout << "active collision: " << meshId << ", "
	//		<< vertexTetMeshId << "\n";

	//	VBDCollisionDetectionResult& collisionResult = getCollisionDetectionResultFromTetMeshId(meshId, vertexTetMeshId);
	//	int iIntersection = 0;
	//	int meshId_intersecting = collisionResult.intersectedTMeshIds[iIntersection];
	//	VBDBaseTetMesh* pTetMesh = tMeshes[meshId_intersecting].get();


	//	VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();
	//	const int surfaceFaceId = collisionResult.closestSurfaceFaceId[iIntersection];

	//	std::cout << "from parallel group: " << pTetMesh->vertexParallelGroups(vertexTetMeshId);
	//	for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
	//	{
	//		IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);
	//		std::cout << ", " << pTetMesh_intersecting->vertexParallelGroups(faceVId);
	//	}
	//	std::cout << "\n";
	//}

	//for (size_t iGroup = 0; iGroup < activeColllisionList.numActiveCollisionsEachParallelGroup.size(); iGroup++)
	//{
	//	size_t numActiveCollision = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
	//	std::cout << "numActiveCollision in group " << iGroup << ": " << numActiveCollision << "\n";
	//	for (size_t iActiveCol = 0; iActiveCol < numActiveCollision; iActiveCol++)
	//	{
	//		int iMesh = activeColllisionList.activeCollisionsEachParallelGroup[iGroup](iActiveCol * 2);
	//		int vId = activeColllisionList.activeCollisionsEachParallelGroup[iGroup](iActiveCol * 2 + 1);
	//		std::cout << "active collision: " << iMesh << ", " << vId << "\n";
	//	}
	//}
	//std::cout << "numActiveDCD: "  << numActiveDCD << "\n";
	//std::cout << "numActiveCCD: "  << numActiveCCD << "\n";
}

void GAIA::VBDPhysics::prepareCollisionDataGPU()
{
	// clear the collision relations
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		for (size_t iVert = 0; iVert < pTetMesh->numVertices(); iVert++)
		{
			CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(iVert);
			collisionDataCPUBuffer.numCollisionRelations = 0;
			collisionDataCPUBuffer.activeColliding = 0;
		}
	}
	activeColllisionList.clear();

	// record the activate collisions
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++) {
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		pTetMesh->activeCollisionMask.setZero();
		std::vector<VBDCollisionDetectionResult>& dcdResults = collisionResultsAll[iMesh];
		if (!pTetMesh->activeForCollision)
		{
			continue;
		}
		// todo: make this parallel using atomic operation
		// loop through all the v-f collisions
		for (int iSurfaceV = 0; iSurfaceV < pTetMesh->surfaceVIds().size(); ++iSurfaceV)
		{
			VBDCollisionDetectionResult& collisionResult = dcdResults[iSurfaceV];
			int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
			CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(surfaceVIdTetMesh);

			if (collisionResult.numIntersections())
			{
				activeColllisionList.addToActiveCollisionList(collisionResult, 1);
				int iIntersection = 0;

				const CollidingPointInfo& collidingPt = collisionResult.collidingPts[iIntersection];
				if (collidingPt.shortestPathFound)
				{
					pTetMesh->activeCollisionMask(surfaceVIdTetMesh) = true;

					int meshId_intersecting = collidingPt.intersectedMeshId;
					int surfaceFaceId = collidingPt.closestSurfaceFaceId;

					collisionDataCPUBuffer.activeColliding = 1;

					collisionDataCPUBuffer.closestSurfaceFaceId = surfaceFaceId;
					collisionDataCPUBuffer.intersectedTMeshId = meshId_intersecting;

					VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();

					for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
					{
						IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);

						collisionDataCPUBuffer.closestSurfaceFaceVIds[iSurfaceFaceV] = faceVId;

						int numColRelations = pTetMesh_intersecting->pTetMeshSharedBase->numCollisionRelation(faceVId);

						if (numColRelations < GPU_COLLISION_RELATION_PREALLOCATION_SIZE)
						{

							//std::cout << "f side: collision mask of " << meshId_intersecting << "," << faceVId << " set to true\n";
							CollisionDataGPU& collisionRelationOtherSide =
								pTetMesh_intersecting->pTetMeshSharedBase->getCollisionDataCPUBuffer(faceVId);
							pTetMesh_intersecting->activeCollisionMask(faceVId) = true;
							pTetMesh_intersecting->pTetMeshSharedBase->numCollisionRelation(faceVId)++;
							CollisionRelationGPU& collisionRelation =
								pTetMesh_intersecting->pTetMeshSharedBase->getCollisionRelationCPUBuffer(faceVId, numColRelations);
							collisionRelation.collisionPrimitiveId = surfaceVIdTetMesh;
							collisionRelation.collisionType = 0;
							collisionRelation.collisionVertexOrder = iSurfaceFaceV;
							collisionRelation.meshId = iMesh;
						}
					}
				}
			}
		}// surface verts

		//for (size_t iVert = 0; iVert < pTetMesh->numVertices(); iVert++)
		//{
		//	CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(iVert);
		//	if (collisionDataCPUBuffer.activeColliding)
		//	{
		//		printf("Vert %d of mesh %d collide with mesh: %d, face: %d, %d, %d\n", iVert, iMesh, collisionDataCPUBuffer.intersectedTMeshId,
		//			collisionDataCPUBuffer.closestSurfaceFaceVIds[0],
		//			collisionDataCPUBuffer.closestSurfaceFaceVIds[1],
		//			collisionDataCPUBuffer.closestSurfaceFaceVIds[2]);
		//	}
		//	/*int numCollisionRelations = collisionDataCPUBuffer.numCollisionRelations;
		//	if (numCollisionRelations)
		//	{
		//		printf("Vertex %d has Number Collision Relations: %d\n", iVert, numCollisionRelations);
		//		for (size_t iRelation = 0; iRelation < numCollisionRelations; iRelation++)
		//		{
		//			CollisionRelationGPU& collisionRelation =
		//				pTetMesh->pTetMeshSharedBase->getCollisionRelationCPUBuffer(iVert, iRelation);
		//			printf("collisionPrimitiveId: %d, collisionVertexOrder: %d \n", collisionRelation.collisionPrimitiveId,
		//				collisionRelation.collisionVertexOrder);
		//		}
		//	}*/
		//}

	}// meshes

	//std::cout << "-------------------------------------\n";

	//for (size_t iCol = 0; iCol < activeColllisionList.activeCollisions.size(); iCol++)
	//{
	//	int meshId = activeColllisionList.activeCollisions[iCol].first;
	//	int vertexTetMeshId = activeColllisionList.activeCollisions[iCol].second;
	//	std::cout << "active collision: " << meshId << ", "
	//		<< vertexTetMeshId << "\n";

	//	VBDCollisionDetectionResult& collisionResult = getCollisionDetectionResultFromTetMeshId(meshId, vertexTetMeshId);
	//	int iIntersection = 0;
	//	int meshId_intersecting = collisionResult.intersectedTMeshIds[iIntersection];
	//	VBDBaseTetMesh* pTetMesh = tMeshes[meshId_intersecting].get();


	//	VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();
	//	const int surfaceFaceId = collisionResult.closestSurfaceFaceId[iIntersection];

	//	std::cout << "from parallel group: " << pTetMesh->vertexParallelGroups(vertexTetMeshId);
	//	for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
	//	{
	//		IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);
	//		std::cout << ", " << pTetMesh_intersecting->vertexParallelGroups(faceVId);
	//	}
	//	std::cout << "\n";
	//}

	//for (size_t iGroup = 0; iGroup < activeColllisionList.numActiveCollisionsEachParallelGroup.size(); iGroup++)
	//{
	//	size_t numActiveCollision = activeColllisionList.numActiveCollisionsEachParallelGroup[iGroup];
	//	std::cout << "numActiveCollision in group " << iGroup << ": " << numActiveCollision << "\n";
	//	for (size_t iActiveCol = 0; iActiveCol < numActiveCollision; iActiveCol++)
	//	{
	//		int iMesh = activeColllisionList.activeCollisionsEachParallelGroup[iGroup](iActiveCol * 2);
	//		int vId = activeColllisionList.activeCollisionsEachParallelGroup[iGroup](iActiveCol * 2 + 1);
	//		std::cout << "active collision: " << iMesh << ", " << vId << "\n";
	//	}
	//}
}

void GAIA::VBDPhysics::dcd()
{
	TICK(timeCsmpUpdatingCollisionInfoDCD);
	if (collisionParams().allowDCD)
	{
		bool rebuildDCDTetSceneBVH = (substep == 0) && !(frameId % physicsParams().dcdTetMeshSceneBVHRebuildSteps) && (frameId);
		bool rebuildDCDSurfaceSceneBVH = (substep == 0) && !(frameId % physicsParams().dcdSurfaceSceneBVHRebuildSteps) && (frameId);
		TICK(timeCsmpUpdatingBVHDCD);
		updateDCDBVH(rebuildDCDTetSceneBVH, rebuildDCDSurfaceSceneBVH);
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingBVHDCD);

		// DCD
		TICK(timeCsmpColDetectDCD);
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
			std::vector<VBDCollisionDetectionResult>& collisionResults = collisionResultsAll[iMesh];
			if (!pTetMesh->activeForCollision)
			{
				continue;
			}
			pTetMesh->penetratedMask.setZero();

			cpu_parallel_for(0, pTetMesh->surfaceVIds().size(), [&](int iSurfaceV)
				{
					int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
					if (pTetMesh->vertex(surfaceVIdTetMesh)[GRAVITY_AXIS] < physicsParams().collisionOffHeight
						// && !pTetMesh->penetratedMask(surfaceVIdTetMesh))
						)
					{
						VBDCollisionDetectionResult& colResult = collisionResults[iSurfaceV];
						pDCD->vertexCollisionDetection(surfaceVIdTetMesh, iMesh, &colResult);
						if (colResult.numIntersections())
						{
							pTetMesh->penetratedMask(surfaceVIdTetMesh) = true;
							ClosestPointQueryResult closestPtResult;
							pDCD->closestPointQuery(&colResult, &closestPtResult);

							// int meshId_intersecting = colResult.intersectedTMeshIds[0];
							// TetMeshFEM* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();
							// int surfaceFaceId = colResult.closestSurfaceFaceId[0];
							//std::cout << colResult.numIntersections() << " collision detected for vertex "
							//	<< surfaceVIdTetMesh << " (iSurfaceV : " << iSurfaceV<< " ) "
							//	<< "pTetMesh->tetVertIndicesToSurfaceVertIndices()[surfaceVIdTetMesh]: " << pTetMesh->tetVertIndicesToSurfaceVertIndices()(surfaceVIdTetMesh)
							//	<< " between mesh " << iMesh << " and " << colResult.intersectedTMeshIds[0]
							//	<< " with face " << surfaceFaceId <<
							//	" [" << pTetMesh_intersecting->surfaceFacesTetMeshVIds()(0, surfaceFaceId) << ", "
							//	<< pTetMesh_intersecting->surfaceFacesTetMeshVIds()(1, surfaceFaceId) << ", "
							//	<< pTetMesh_intersecting->surfaceFacesTetMeshVIds()(2, surfaceFaceId)
							//	<< "]\n";
						}
					}

				});
		}
		TOCK_STRUCT(timeStatistics(), timeCsmpColDetectDCD);
	}
	else if (collisionParams().allowCCD)
		// if do ccd only then clear all the collision results
	{
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			TetMeshFEM* pTetMesh = tMeshes[iMesh].get();
			std::vector<VBDCollisionDetectionResult>& collisionResults = collisionResultsAll[iMesh];

			TICK(timeCsmpColDetectDCD);
			cpu_parallel_for(0, pTetMesh->surfaceVIds().size(), [&](int iSurfaceV)
				{
					VBDCollisionDetectionResult& colResult = collisionResults[iSurfaceV];
					colResult.clear();
				});
		}
	}
	TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingCollisionInfoDCD);
}

void GAIA::VBDPhysics::intermediateDCD()
{
	TICK(timeCsmpUpdatingCollisionInfoDCD);
	if (collisionParams().allowDCD)
	{
		bool rebuildDCDTetSceneBVH = false;
		bool rebuildDCDSurfaceSceneBVH = false;
		TICK(timeCsmpUpdatingBVHDCD);
		updateDCDBVH(rebuildDCDTetSceneBVH, rebuildDCDSurfaceSceneBVH);
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingBVHDCD);

		// DCD
		TICK(timeCsmpColDetectDCD);
		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
			std::vector<VBDCollisionDetectionResult>& collisionResults = collisionResultsAll[iMesh];
			if (!pTetMesh->activeForCollision)
			{
				continue;
			}
			cpu_parallel_for(0, pTetMesh->surfaceVIds().size(), [&](int iSurfaceV)
				{
					int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
					if (pTetMesh->penetratedMask(surfaceVIdTetMesh)
						|| !collisionParams().allowCCD)
					{
						VBDCollisionDetectionResult& colResult = collisionResults[iSurfaceV];
						pDCD->vertexCollisionDetection(surfaceVIdTetMesh, iMesh, &colResult);
						if (colResult.numIntersections())
						{
							ClosestPointQueryResult closestPtResult;
							pDCD->closestPointQuery(&colResult, &closestPtResult);

						}
					}
				});
		}
		TOCK_STRUCT(timeStatistics(), timeCsmpColDetectDCD);
	}

	TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingCollisionInfoDCD);
}

void GAIA::VBDPhysics::ccd()
{
	TICK(timeCsmpUpdatingCollisionInfoCCD);
	if (collisionParams().allowCCD)
	{
		TICK(timeCsmpUpdatingBVHCCD);
		bool rebuildCCDBVH = (substep == 0) && !(frameId % physicsParams().ccdBVHRebuildSteps) && (frameId);
		updateCCDBVH(rebuildCCDBVH);
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingBVHCCD);

		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
			std::vector<VBDCollisionDetectionResult>& collisionResults = collisionResultsAll[iMesh];
			if (!pTetMesh->activeForCollision)
			{
				continue;
			}
			TICK(timeCsmpColDetectCCD);
			cpu_parallel_for(0, pTetMesh->surfaceVIds().size(), [&](int iSurfaceV)
				{
					VBDCollisionDetectionResult& colResult = collisionResults[iSurfaceV];
					// if not already penetrated use ccd
					int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
					if (
						pTetMesh->vertex(surfaceVIdTetMesh)[GRAVITY_AXIS] < physicsParams().collisionOffHeight
						&& !pTetMesh->penetratedMask(surfaceVIdTetMesh)
						)
					{
						pCCD->vertexContinuousCollisionDetection(surfaceVIdTetMesh, iMesh, &colResult);

					}
				}
			);
			TOCK_STRUCT(timeStatistics(), timeCsmpColDetectCCD);
		}
	}
	TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingCollisionInfoCCD);
}

void GAIA::VBDPhysics::intermediateCCD()
{
	TICK(timeCsmpUpdatingCollisionInfoCCD);
	if (collisionParams().allowCCD)
	{
		TICK(timeCsmpUpdatingBVHCCD);
		bool rebuildCCDBVH = false;
		updateCCDBVH(rebuildCCDBVH);
		TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingBVHCCD);

		for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
		{
			VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
			std::vector<VBDCollisionDetectionResult>& collisionResults = collisionResultsAll[iMesh];
			if (!pTetMesh->activeForCollision)
			{
				continue;
			}
			TICK(timeCsmpColDetectCCD);
			cpu_parallel_for(0, pTetMesh->surfaceVIds().size(), [&](int iSurfaceV)
				{
					VBDCollisionDetectionResult& colResult = collisionResults[iSurfaceV];
					int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
					// if not already penetrated use ccd
					if (
						pTetMesh->vertex(surfaceVIdTetMesh)[GRAVITY_AXIS] < physicsParams().collisionOffHeight
						&& !pTetMesh->penetratedMask(surfaceVIdTetMesh)
						)
					{
						pCCD->vertexContinuousCollisionDetection(surfaceVIdTetMesh, iMesh, &colResult);

					}
				}
			);
			TOCK_STRUCT(timeStatistics(), timeCsmpColDetectCCD);
		}
	}
	TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingCollisionInfoCCD);
}

void GAIA::VBDPhysics::intermediateCollisionDetection()
{
	if (physicsParams().useGPU) {
		syncAllToCPUVertPosOnly(true);
	}
	intermediateDCD();
	intermediateCCD();
	if (physicsParams().useGPU) {
		prepareCollisionDataGPU();
		syncAllToGPU(false);
	}
	else {
		prepareCollisionDataCPU();
	}
}

void GAIA::VBDPhysics::updateDCDBVH(bool rebuildTetMeshScene, bool rebuildSurfaceScene)
{
	RTCBuildQuality tetSceneQuality = rebuildTetMeshScene ? RTC_BUILD_QUALITY_LOW : RTC_BUILD_QUALITY_REFIT;
	RTCBuildQuality surfaceSceneQuality = rebuildSurfaceScene ? RTC_BUILD_QUALITY_LOW : RTC_BUILD_QUALITY_REFIT;
	pDCD->updateBVH(tetSceneQuality, surfaceSceneQuality, true);
}

void GAIA::VBDPhysics::updateCCDBVH(bool rebuildScene)
{
	TICK(timeCsmpUpdatingBVHCCD);
	RTCBuildQuality sceneQuality = rebuildScene ? RTC_BUILD_QUALITY_LOW : RTC_BUILD_QUALITY_REFIT;
	pCCD->updateBVH(sceneQuality);
	TOCK_STRUCT(timeStatistics(), timeCsmpUpdatingBVHCCD);
}

void GAIA::VBDPhysics::updateAllCollisionInfos()
{
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		cpu_parallel_for(0, pTetMesh->numVertices(), [&](int iV)
			{
				if (pTetMesh->activeCollisionMask[iV])
				{
					VBDCollisionDetectionResult& colResult = getCollisionDetectionResultFromTetMeshId(iMesh, iV);
					updateCollisionInfo(colResult);
				}
			});
	}
}

void GAIA::VBDPhysics::updateCollisionInfo(VBDCollisionDetectionResult& collisionResult)
{
	if (collisionResult.collisionForceAndHessian.size() < collisionResult.numIntersections())
	{
		collisionResult.collisionForceAndHessian.insertN(collisionResult.collisionForceAndHessian.size() - 1, 
			collisionResult.numIntersections() - collisionResult.collisionForceAndHessian.size());
	}
	
	for (size_t iIntersection = 0; iIntersection < collisionResult.numIntersections(); iIntersection++)
	{

		CollidingPointInfo& collidingPt = collisionResult.collidingPts[iIntersection];

		if (collidingPt.shortestPathFound)
		{
			int curMeshID = collisionResult.idTMQuery;
			int collidedMeshID = collidingPt.intersectedMeshId;
			int closestFaceId = collidingPt.closestSurfaceFaceId;

			VBDBaseTetMesh* pCurTM = tMeshes[curMeshID].get();
			VBDBaseTetMesh* pIntersectedTM = tMeshes[collidedMeshID].get();

			int t0 = pIntersectedTM->surfaceFacesTetMeshVIds()(0, closestFaceId);
			int t1 = pIntersectedTM->surfaceFacesTetMeshVIds()(1, closestFaceId);
			int t2 = pIntersectedTM->surfaceFacesTetMeshVIds()(2, closestFaceId);

			embree::Vec3fa a = loadVertexPos(pIntersectedTM, t0),
				b = loadVertexPos(pIntersectedTM, t1),
				c = loadVertexPos(pIntersectedTM, t2),
				p = loadVertexPos(pCurTM, collisionResult.idVQuery);

			embree::Vec3fa bary;
			ClosestPointOnTriangleType closestPtType;
			embree::Vec3fa closestSurfacePoint_ = closestPointTriangle(p, a, b, c, bary, closestPtType);

			collidingPt.closestSurfacePt << closestSurfacePoint_.x, closestSurfacePoint_.y, closestSurfacePoint_.z;
			collidingPt.closestSurfacePtBarycentrics << bary.x, bary.y, bary.z;

			Vec3 normal;
			Vec3 faceNormal;
			embree::Vec3fa closestP_to_p;

			pIntersectedTM->computeFaceNormal(closestFaceId, faceNormal);

			if (closestPtType == GAIA::ClosestPointOnTriangleType::AtInterior)
			{
				collidingPt.closestPointNormal = faceNormal;

			}
			// all other cases
			else if (closestPtType != GAIA::ClosestPointOnTriangleType::NotFound)
			{
				closestP_to_p = closestSurfacePoint_ - p;
				closestP_to_p /= embree::length(closestP_to_p);
				normal << closestP_to_p.x, closestP_to_p.y, closestP_to_p.z;

				if (normal.dot(faceNormal) < 0)
				{
					normal = -normal;;
				}
				collidingPt.closestPointNormal = normal;
			}

			// Frictionand Collision
			// Get the relative displacement

			VBDCollisionInfo& collisionForceAndHessian = collisionResult.collisionForceAndHessian[iIntersection];

			Vec3 x = pCurTM->vertex(collisionResult.idVQuery);
			collisionForceAndHessian.diff = collidingPt.closestSurfacePt - x;
			Vec3 dx3 = x - pCurTM->vertexPrevPos(collisionResult.idVQuery);
			Vec3 dx0 = pIntersectedTM->vertex(t0) - pIntersectedTM->vertexPrevPos(t0);
			Vec3 dx1 = pIntersectedTM->vertex(t1) - pIntersectedTM->vertexPrevPos(t1);
			Vec3 dx2 = pIntersectedTM->vertex(t2) - pIntersectedTM->vertexPrevPos(t2);
			Vec3 dx = dx3 - (bary.x * dx0 + bary.y * dx1 + bary.z * dx2);
			Vec3 e0 = (pIntersectedTM->vertex(t1) - pIntersectedTM->vertex(t0)).normalized();
			Mat3x2 T = Mat3x2::Zero();
			T.col(0) = e0;
			const Vec3& n = collidingPt.closestPointNormal;
			T.col(1) = (e0.cross(n)).normalized();
			Vec2 u = T.transpose() * dx;

			// average of the two friction coefficients
			CFloatingType mu = (getObjectParam(collidedMeshID).frictionDynamic + getObjectParam(curMeshID).frictionDynamic) * 0.5;
			CFloatingType epsV = (getObjectParam(collidedMeshID).frictionEpsV + getObjectParam(curMeshID).frictionEpsV) * 0.5;
			CFloatingType dt = physicsParams().dt;
			CFloatingType epsU = epsV * dt;
			//Vec3 frictionForce;
			//Mat3 frictionForceHessian;
			computeVertexFriction(mu, 1.0f, T, u, epsU, collisionForceAndHessian.frictionForce, collisionForceAndHessian.frictionHessian);

		}
	}
}

void GAIA::VBDPhysics::solveCollisionsSequentially()
{
	for (size_t iMesh = 0; iMesh < basetetMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		for (int iV = 0; iV < pTetMesh->numVertices(); ++iV)
		{
			if (pTetMesh->activeCollisionMask[iV] && pTetMesh->activeForMaterialSolve)
			{
				updateCollisionInfo(getCollisionDetectionResultFromTetMeshId(iMesh, iV));
				VBDStepWithCollision(pTetMesh, iMesh, iV);
			}
		}
	}

}

void GAIA::VBDPhysics::VBDStep(TetMeshFEM* pMesh_, IdType vertexId)
{

	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)pMesh_;

	Mat3 h;
	Vec3 force;
	h.setZero();
	force.setZero();

	pMesh->accumlateInertiaForceAndHessian(vertexId, force, h);
	//printf("-------------------------------------\nProcessing vertex %d\n", vertexId);
	//printf("after accumulateInertiaForceAndHessian\n");
	//printf("force: ");
	//CuMatrix::printFloatVec(force.data(), 3);
	//printf("h: ");
	//CuMatrix::printMat3(h.data());
	Mat3 tmp_h = h;
	pMesh->accumlateMaterialForceAndHessian2(vertexId, force, h);
	Mat3 K = h - tmp_h;
	// CuMatrix::printMat3(K.data());
	accumlateDampingForceAndHessian(pMesh, vertexId, force, h, K);
	//printf("after accumlateMaterialForceAndHessian\n");
	//printf("force: ");
	//CuMatrix::printFloatVec(force.data(), 3);
	//printf("h: ");
	//CuMatrix::printMat3(h.data());

	accumlateBoundaryForceAndHessian(pMesh, 0, vertexId, force, h);
	//accum
	// add external force
	force += pMesh->vertexExternalForces.col(vertexId);

	// line search
	if (force.squaredNorm() > CMP_EPSILON2)
	{
		// Solve linear system
		Vec3 descentDirection;
		FloatingType stepSize = physicsParams().stepSize;
		FloatingType lineSearchShrinkFactor = physicsParams().tau;
		bool solverSuccess = CuMatrix::solve3x3_psd_stable(h.data(), force.data(), descentDirection.data());

		if (!solverSuccess)
		{
			stepSize = physicsParams().stepSizeGD;
			descentDirection = force;
			lineSearchShrinkFactor = physicsParams().tau_GD;
			std::cout << "Solver failed at vertex " << vertexId << std::endl;
		}

		// line search
#ifdef APPLY_LOCAL_LINE_SEARCH
		FloatingType meInertia = 0;
		FloatingType meElastic_elastic = 0;
		FloatingType initialEnergy = pMesh->evaluateVertexMeritEnergy(vertexId, meInertia, meElastic_elastic);
		FloatingType e = pMesh->backTracingLineSearchVBD(vertexId, descentDirection, initialEnergy, stepSize, 0.f,
			lineSearchShrinkFactor, physicsParams().backtracingLineSearchMaxIters);

		if (isnan(e))
		{
			assert(false);
		}
#else
		//printf("descentDirection: ");
		//CuMatrix::printFloatVec(descentDirection.data(), 3);

		pMesh->vertex(vertexId) += stepSize * descentDirection;
#endif // LOCAL_LINE_SEARCH

	}
}

void GAIA::VBDPhysics::VBDStepWithCollision(TetMeshFEM* pMesh_, IdType meshId, IdType vertexId, bool apply_friction)
{
	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)pMesh_;

	Mat3 h;
	Vec3 force;
	h.setZero();
	force.setZero();

	pMesh->accumlateMaterialForceAndHessian2(vertexId, force, h);
	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
		MaterialForce[meshId].col(vertexId) += force;
		});
	pMesh->accumlateInertiaForceAndHessian(vertexId, force, h);
	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
		InertiaForce[meshId].col(vertexId) += force - MaterialForce[meshId].col(vertexId);
		});
	/*Mat3 tmp_h = h;
	pMesh->accumlateMaterialForceAndHessian(vertexId, force, h);
	Mat3 K = h - tmp_h;
	accumlateDampingForceAndHessian(pMesh, vertexId, force, h, K);*/
	accumlateBoundaryForceAndHessian(pMesh, meshId, vertexId, force, h, apply_friction);
	accumlateCollisionForceAndHessian(pMesh, meshId, vertexId, force, h, apply_friction);

	//accum
	// add external force
	force += pMesh->vertexExternalForces.col(vertexId);

	// line search
	if (force.squaredNorm() > CMP_EPSILON2)
	{
		// Solve linear system
		Vec3 descentDirection;
		FloatingType stepSize = physicsParams().stepSize;
		FloatingType lineSearchShrinkFactor = physicsParams().tau;

		bool solverSuccess;
		if (physicsParams().useDouble3x3) {
			//construct a double version of h, forcem and positionsNew
			double H[9] = { h(0,0), h(1,0), h(2,0),
				h(0,1), h(1,1), h(2,1), h(0,2), h(1,2), h(2,2) };
			double F[3] = { force(0), force(1), force(2) };
			double dx[3] = { 0, 0, 0 };
			solverSuccess = CuMatrix::solve3x3_psd_stable(H, F, dx);
			descentDirection = Vec3(dx[0], dx[1], dx[2]);
		}
		else {
			solverSuccess = CuMatrix::solve3x3_psd_stable(h.data(), force.data(), descentDirection.data());
		}

		if (!solverSuccess)
		{
			stepSize = physicsParams().stepSizeGD;
			descentDirection = force;
			lineSearchShrinkFactor = physicsParams().tau_GD;
			std::cout << "Solver failed at vertex " << vertexId << std::endl;
			// std::exit(-1);
		}
		if (descentDirection.hasNaN()) {
			std::cout << "force: " << force.transpose() << "\nHessian:\n" << h;
			std::cout << "descentDirection has NaN at vertex " << vertexId << std::endl;
			std::exit(-1);
		}

		// line search
#ifdef APPLY_LOCAL_LINE_SEARCH
		FloatingType meInertia = 0;
		FloatingType meElastic_elastic = 0;
		FloatingType initialEnergy = pMesh->evaluateVertexMeritEnergy(vertexId, meInertia, meElastic_elastic);
		FloatingType e = pMesh->backTracingLineSearchVBD(vertexId, descentDirection, initialEnergy, stepSize, 0.f,
			lineSearchShrinkFactor, physicsParams().backtracingLineSearchMaxIters);

		if (isnan(e))
		{
			assert(false);
		}
#else

		pMesh->vertex(vertexId) += stepSize * descentDirection;
#endif // LOCAL_LINE_SEARCH
	}
}


void GAIA::VBDPhysics::updateVelocities()
{
	// updateAllCollisionInfos();
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDBaseTetMesh* pMesh = tMeshes[iMesh].get();

		if (pMesh->activeForMaterialSolve)
		{
			cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {
				updateVelocity(pMesh, iV);
				// applyBoudnaryFriction(pMesh, iV);
				// applyCollisionFriction(pMesh, iMesh, iV);
				});
		}

	}

	dampVelocities();

}

void GAIA::VBDPhysics::updateVelocitiesGPU()
{
	VBDUpdateVelocityGPU(getVBDPhysicsDataGPU(), vertexAllParallelGroupsBuffer->getGPUBuffer(), numAllVertices,
		physicsParams().numThreadsVBDSolve, cudaStream);
}

void GAIA::VBDPhysics::updateVelocity(TetMeshFEM* pMesh, IdType vertexId)
{
	pMesh->velocity(vertexId) = (pMesh->vertex(vertexId) - pMesh->vertexPrevPos(vertexId)) / physicsParams().dt;

}

void GAIA::VBDPhysics::dampVelocities()
{
	for (size_t iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		if (!tMeshes[iMesh]->activeForMaterialSolve)
		{
			continue;
		}
		//for (size_t iVert = 0; iVert < tMeshes[iMesh]->numVertices(); iVert++)
		cpu_parallel_for(0, tMeshes[iMesh]->numVertices(), [&](int iVert)
			{
				FloatingType vMag = tMeshes[iMesh]->mVelocity.col(iVert).norm();
				FloatingType maxVelocityMagnitude = getObjectParam(iMesh).maxVelocityMagnitude;

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

		for (size_t iFixedPoint = 0; iFixedPoint < tMeshes[iMesh]->pObjectParams->fixedPoints.size(); iFixedPoint++)
		{
			tMeshes[iMesh]->velocity(tMeshes[iMesh]->pObjectParams->fixedPoints[iFixedPoint]).setZero();
		}
	}
}

void GAIA::VBDPhysics::computeVertexFriction(CFloatingType mu, CFloatingType lambda, const Mat3x2& T, const Vec2& u, CFloatingType epsU, Vec3& force, Mat3& hessian)
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
		//hessian = mu * lambda * T * (df1_x_minus_f1_over_x3 * u * u.transpose() + f1_SF_over_x * Mat2::Identity()) * T.transpose();
		hessian = mu * lambda * T * (f1_SF_over_x * Mat2::Identity()) * T.transpose();
	}
	else
	{
		force.setZero();
		hessian.setZero();
	}
}

void GAIA::VBDPhysics::applyBoudnaryFriction(TetMeshFEM* pMesh, IdType vertexId)
{
	CFloatingType boundaryCollisionStiffness = physicsParams().boundaryCollisionStiffness;
	CFloatingType dt = physicsParams().dt;
	CFloatingType vertInvMass = pMesh->vertexInvMass(vertexId);
	CFloatingType frictionRatio = physicsParams().boundaryFrictionDynamic;
	for (size_t iDim = 0; iDim < 3; iDim++)
	{
		Vec3 contactNormal = Vec3::Zero();

		CFloatingType lowerBound = physicsParams().worldBounds(iDim, 0);
		CFloatingType upperBound = physicsParams().worldBounds(iDim, 1);

		if (pMesh->vertex(vertexId)[iDim] < lowerBound)
		{
			CFloatingType penetrationDepth = lowerBound - pMesh->vertex(vertexId)[iDim];

			contactNormal(iDim) = 1;
			pMesh->velocity(vertexId) =
				applyFrictinalVelocityDamping(pMesh->velocity(vertexId), contactNormal, vertInvMass, frictionRatio, penetrationDepth * boundaryCollisionStiffness, dt);
		}
		else if (pMesh->vertex(vertexId)[iDim] > upperBound)
		{
			CFloatingType penetrationDepth = pMesh->vertex(vertexId)[iDim] - upperBound;

			contactNormal(iDim) = -1;
			pMesh->velocity(vertexId) =
				applyFrictinalVelocityDamping(pMesh->velocity(vertexId), contactNormal, vertInvMass, frictionRatio, penetrationDepth * boundaryCollisionStiffness, dt);
		}
	}
}

void GAIA::VBDPhysics::applyCollisionFriction(VBDBaseTetMesh* pMesh, IdType meshId, IdType vertexId)
{
	if (pMesh->activeCollisionMask[vertexId])
	{
		CFloatingType dt = physicsParams().dt;
		CFloatingType vertInvMass = pMesh->vertexInvMass(vertexId);
		CFloatingType frictionRatio = pMesh->pObjParamsVBD->frictionDynamic;

		int surfaceVId = pMesh->tetVertIndicesToSurfaceVertIndices()(vertexId);
		FloatingType collisionForce = 0.f;
		VBDCollisionDetectionResult& colResult = getCollisionDetectionResultFromSurfaceMeshId(meshId, surfaceVId);
		Vec3 contactNormal;
		// v side of the v-f collisions
		if (colResult.numIntersections())
		{
			for (size_t iIntersection = 0; iIntersection < colResult.numIntersections(); iIntersection++)
			{
				collisionForce = computeCollisionRepulsiveForcePerCollision(iIntersection, 3, colResult, contactNormal);
				pMesh->velocity(vertexId) =
					applyFrictinalVelocityDamping(pMesh->velocity(vertexId), contactNormal, vertInvMass, frictionRatio, collisionForce, dt);
			}
		}

		// f side of the v-f collision
		CollisionRelationList& colRelations = getCollisionRelationList(meshId, vertexId);
		for (size_t iCollision = 0; iCollision < colRelations.size(); iCollision++)
		{
			CollisionRelation& colRelation = colRelations[iCollision];

			VBDCollisionDetectionResult& colResultOther = getCollisionDetectionResultFromSurfaceMeshId(colRelation.meshId, colRelation.surfaceVertexId);
			collisionForce = computeCollisionRepulsiveForcePerCollision(colRelation.collisionId, colRelation.collisionVertexOrder, colResultOther, contactNormal);
			pMesh->velocity(vertexId) =
				applyFrictinalVelocityDamping(pMesh->velocity(vertexId), contactNormal, vertInvMass, frictionRatio, collisionForce, dt);
		}
	}

}

//FloatingType GAIA::VBDPhysics::computeCollisionRepulsiveForce(VBDBaseTetMesh* pMesh, IdType meshId, IdType vertexId)
//{
//	int surfaceVId = pMesh->tetVertIndicesToSurfaceVertIndices()(vertexId);
//	FloatingType collisionForce = 0.f;
//	VBDCollisionDetectionResult& colResult = getCollisionDetectionResult(meshId, surfaceVId);
//	// v side of the v-f collisions
//	if (colResult.numIntersections())
//	{
//		for (size_t iIntersection = 0; iIntersection < colResult.numIntersections(); iIntersection++)
//		{
//			collisionForce += computeCollisionRepulsiveForcePerCollision(iIntersection, 3, colResult);
//		}
//	}
//
//	// f side of the v-f collision
//	CollisionRelationList& colRelations = getCollisionRelationList(meshId, vertexId);
//	for (size_t iCollision = 0; iCollision < colRelations.size(); iCollision++)
//	{
//		CollisionRelation& colRelation = colRelations[iCollision];
//
//		VBDCollisionDetectionResult& colResultOther = getCollisionDetectionResult(colRelation.meshId, colRelation.surfaceVertexId);
//		collisionForce += computeCollisionRepulsiveForcePerCollision(colRelation.collisionId, colRelation.collisionVertexOrder, colResultOther);
//
//	}
//	return collisionForce;
//}

FloatingType GAIA::VBDPhysics::computeCollisionRepulsiveForcePerCollision(int collisionId, int collisionVertexOrder, VBDCollisionDetectionResult& collisionResult,
	Vec3& contactNormal)
{
	int curMeshID = collisionResult.idTMQuery;
	VBDBaseTetMesh* pCurTM = tMeshes[curMeshID].get();
	CollidingPointInfo& collidingPt = collisionResult.collidingPts[collisionId];

	if (!collidingPt.shortestPathFound)
	{
		return 0.f;
	}

	int collidedMeshID = collidingPt.intersectedMeshId;
	int closestFaceId = collidingPt.closestSurfaceFaceId;
	VBDBaseTetMesh* pIntersectedTM = tMeshes[collidedMeshID].get();

	FloatingType b;
	Vec3 p = pCurTM->vertex(collisionResult.idVQuery);
	Vec3 diff, c, n;
	CFloatingType k = physicsParams().collisionStiffness;
	CFloatingType collisionRadius = physicsParams().collisionAirDistance;
	n = collidingPt.closestPointNormal;
	c = collidingPt.closestSurfacePt;
	diff = p - c;

	contactNormal = n;
	FloatingType penetrationDepth = -diff.dot(n);

	if (physicsParams().collisionEnergyType == 1)
	{
		penetrationDepth += physicsParams().collisionAirDistance;
	}

	FloatingType collisionRepulsiveForce = 0.f;
	if (penetrationDepth > 0)
	{
		switch (collisionVertexOrder)
		{
		case 0:
			b = -collidingPt.closestSurfacePtBarycentrics[0];
			break;
		case 1:
			b = -collidingPt.closestSurfacePtBarycentrics[1];
			break;
		case 2:
			b = -collidingPt.closestSurfacePtBarycentrics[2];
			break;
		case 3:
			b = 1.0f;
			// diff = -diff;
			break;
		default:
			break;
		}

		if (b > 0)
		{
			switch (physicsParams().collisionEnergyType)
			{
			case 0:
				collisionRepulsiveForce = k * b * diff.norm();
				break;
			case 1:
				// simplified point-plane
				collisionRepulsiveForce = k * b * penetrationDepth;
				break;
			default:
				break;
			}

		}
	}
	return abs(collisionRepulsiveForce);
}

Vec3 GAIA::VBDPhysics::applyFrictinalVelocityDamping(Vec3 velocity, Vec3& contactNormal, CFloatingType vertInvMass,
	CFloatingType frictionRatio, CFloatingType contactForce, CFloatingType dt)
{
	Vec3 contactVel = contactNormal * (velocity.dot(contactNormal));
	Vec3 orthogonalVel = velocity - contactVel;
	FloatingType orthogonalVelNorm = orthogonalVel.norm();
	FloatingType velDamping = contactForce * vertInvMass * dt * frictionRatio;

	if (velDamping >= orthogonalVelNorm)
	{
		orthogonalVel.setZero();
	}
	else {
		FloatingType orthogonalVelNormDamped = orthogonalVelNorm - velDamping;
		orthogonalVel = orthogonalVel * (orthogonalVelNormDamped / orthogonalVelNorm);
	}

	return orthogonalVel + contactVel;
}

void GAIA::VBDPhysics::accumlateDampingForceAndHessian(TetMeshFEM* pMesh, IdType vertexId, Vec3& force, Mat3& hessian, const Mat3& K)
{
	CFloatingType dt = physicsParams().dt;
	// gamma is the damping ratio
	CFloatingType gamma = pMesh->pObjectParams->dampingGamma;

	Vec3 velocity = (pMesh->vertex(vertexId) - pMesh->vertexPrevPos(vertexId)) / dt;
	Vec3 dampingForce = -gamma * (K * velocity);
	Mat3 dampingHessian = (gamma / dt) * K;

	force += dampingForce;
	hessian += dampingHessian;
}


void GAIA::VBDPhysics::solveBoxBoundaryConstraintForVertex(TetMeshFEM* pMesh, IdType vertexId)
{
	for (size_t iDim = 0; iDim < 3; iDim++)
	{
		Vec3 contactNormal = Vec3::Zero();

		CFloatingType lowerBound = physicsParams().worldBounds(iDim, 0);
		CFloatingType upperBound = physicsParams().worldBounds(iDim, 1);

		if (pMesh->vertex(vertexId)[iDim] < lowerBound)
		{
			pMesh->vertex(vertexId)[iDim] = lowerBound;
		}
		else if (pMesh->vertex(vertexId)[iDim] > upperBound)
		{
			pMesh->vertex(vertexId)[iDim] = upperBound;
		}
	}
}

void GAIA::VBDPhysics::accumlateBoundaryForceAndHessian(TetMeshFEM* pMesh, IdType meshId, IdType vertexId, Vec3& force, Mat3& hessian, bool apply_friction)
{
	FloatingType boundaryCollisionStiffness = physicsParams().boundaryCollisionStiffness;
	if (physicsParams().useBowlGround) {
		const Vec3& center = physicsParams().bowlCenter;
		const FloatingType& radius = physicsParams().bowlRadius;
		const Vec3 pos = pMesh->vertex(vertexId);
		// If no cap and upper half plane, ignore
		if (!physicsParams().bowlCap && (pos.y() - center.y() > 0))return;
		// If inside the sphere, ignore
		CFloatingType dist = (pos - center).norm();
		if (dist < radius)return;
		const Vec3 dir = (pos - center) / dist;

		// Penalty Force
		CFloatingType penetrationDepth = dist - radius;
		debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
			boundaryCollisionForce[meshId].col(vertexId) += -penetrationDepth * boundaryCollisionStiffness * dir;
			});
		force += -penetrationDepth * boundaryCollisionStiffness * dir;
		hessian += boundaryCollisionStiffness * (-Mat3::Identity() + dir * dir.transpose() + radius / dist * (Mat3::Identity() - dir * dir.transpose()));

		// Friction
		if (apply_friction) {
			Vec3 dx = pMesh->vertex(vertexId) - pMesh->vertexPrevPos(vertexId);
			CFloatingType dt = physicsParams().dt;
			Mat3x2 T = Mat3x2::Zero();
			T.col(0) = (dir.cross(Vec3::UnitZ())).normalized();
			T.col(1) = (dir.cross(T.col(0))).normalized();
			Vec2 u = T.transpose() * dx;
			CFloatingType lambda = penetrationDepth * boundaryCollisionStiffness;
			CFloatingType mu = physicsParams().boundaryFrictionDynamic;
			CFloatingType epsV = physicsParams().boundaryFrictionEpsV;
			CFloatingType epsU = epsV * dt;
			Vec3 frictionForce;
			Mat3 frictionForceHessian;
			computeVertexFriction(mu, lambda, T, u, epsU, frictionForce, frictionForceHessian);
			debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
				boundaryFrictionForce[meshId].col(vertexId) += frictionForce;
				});
			force += frictionForce;
			hessian += frictionForceHessian;
		}
	}
	else if (physicsParams().usePlaneGround) {
		for (size_t iDim = 0; iDim < 3; iDim++)
		{
			Vec3 contactNormal = Vec3::Zero();
			CFloatingType dt = physicsParams().dt;

			CFloatingType lowerBound = physicsParams().worldBounds(iDim, 0);
			CFloatingType upperBound = physicsParams().worldBounds(iDim, 1);

			if (pMesh->vertex(vertexId)[iDim] < lowerBound)
			{
				// Penalty Force
				CFloatingType penetrationDepth = lowerBound - pMesh->vertex(vertexId)[iDim];
				debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
					boundaryCollisionForce[meshId](iDim, vertexId) += penetrationDepth * boundaryCollisionStiffness;
					});
				force(iDim) += penetrationDepth * boundaryCollisionStiffness;
				hessian(iDim, iDim) += boundaryCollisionStiffness;

				// Friction
				if (apply_friction) {
					Vec3 dx = pMesh->vertex(vertexId) - pMesh->vertexPrevPos(vertexId);
					Mat3x2 T = Mat3x2::Zero();
					T((iDim + 1) % 3, 0) = 1;
					T((iDim + 2) % 3, 1) = 1;
					Vec2 u = T.transpose() * dx;
					CFloatingType lambda = penetrationDepth * boundaryCollisionStiffness;
					CFloatingType mu = physicsParams().boundaryFrictionDynamic;
					CFloatingType epsV = physicsParams().boundaryFrictionEpsV;
					CFloatingType epsU = epsV * dt;
					Vec3 frictionForce;
					Mat3 frictionForceHessian;
					computeVertexFriction(mu, lambda, T, u, epsU, frictionForce, frictionForceHessian);
					//CFloatingType frictionForceNorm = frictionForce.norm();
					//CFloatingType maxFrictionForceNorm = u.norm() * pMesh->vertexMass(vertexId);
					//if (frictionForceNorm > maxFrictionForceNorm)
					//{
					//	CFloatingType ratio = maxFrictionForceNorm / frictionForceNorm;
					//	// frictionForce *= ratio;
					//	// frictionForceHessian *= ratio;
					//}
					debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
						boundaryFrictionForce[meshId].col(vertexId) += frictionForce;
						});
					force += frictionForce;
					hessian += frictionForceHessian;
				}

			}
			else if (pMesh->vertex(vertexId)[iDim] > upperBound)
			{
				CFloatingType penetrationDepth = pMesh->vertex(vertexId)[iDim] - upperBound;
				debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
					boundaryCollisionForce[meshId](iDim, vertexId) -= penetrationDepth * boundaryCollisionStiffness;
					});
				force(iDim) -= penetrationDepth * boundaryCollisionStiffness;
				hessian(iDim, iDim) += boundaryCollisionStiffness;

				// Friction
				if (apply_friction) {
					Vec3 dx = pMesh->vertex(vertexId) - pMesh->vertexPrevPos(vertexId);
					Mat3x2 T = Mat3x2::Zero();
					T((iDim + 1) % 3, 0) = 1;
					T((iDim + 2) % 3, 1) = 1;
					Vec2 u = T.transpose() * dx;
					CFloatingType lambda = penetrationDepth * boundaryCollisionStiffness;
					CFloatingType mu = physicsParams().boundaryFrictionDynamic;
					CFloatingType epsV = physicsParams().boundaryFrictionEpsV;
					CFloatingType epsU = epsV * dt;
					Vec3 frictionForce;
					Mat3 frictionForceHessian;
					computeVertexFriction(mu, lambda, T, u, epsU, frictionForce, frictionForceHessian);
					//CFloatingType frictionForceNorm = frictionForce.norm();
					//CFloatingType maxFrictionForceNorm = u.norm() * pMesh->vertexMass(vertexId);
					//if (frictionForceNorm > maxFrictionForceNorm)
					//{
					//	CFloatingType ratio = maxFrictionForceNorm / frictionForceNorm;
					//	// frictionForce *= ratio;
					//	// frictionForceHessian *= ratio;
					//}
					debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
						boundaryFrictionForce[meshId].col(vertexId) += frictionForce;
						});
					force += frictionForce;
					hessian += frictionForceHessian;
				}
			}
		}
	}
}

void GAIA::VBDPhysics::accumlateCollisionForceAndHessian(TetMeshFEM* pMesh, IdType meshId, IdType vertexId, Vec3& force, Mat3& hessian, bool apply_friction)
{
	int surfaceVId = pMesh->tetVertIndicesToSurfaceVertIndices()(vertexId);
	if (surfaceVId >= 0)
	{
		VBDCollisionDetectionResult& colResult = getCollisionDetectionResultFromSurfaceMeshId(meshId, surfaceVId);
		// v side of the v-f collisions
		if (colResult.numIntersections())
		{
			for (size_t iIntersection = 0; iIntersection < colResult.numIntersections(); iIntersection++)
			{
				accumlateCollisionForceAndHessianPerCollision(iIntersection, 3, colResult, force, hessian, meshId, vertexId, apply_friction);
			}
		}

		// f side of the v-f collision
		CollisionRelationList& colRelations = getCollisionRelationList(meshId, vertexId);
		for (size_t iCollision = 0; iCollision < colRelations.size(); iCollision++)
		{
			CollisionRelation& colRelation = colRelations[iCollision];

			VBDCollisionDetectionResult& colResultOther = getCollisionDetectionResultFromSurfaceMeshId(colRelation.meshId, colRelation.surfaceVertexId);
			accumlateCollisionForceAndHessianPerCollision(colRelation.collisionId, colRelation.collisionVertexOrder, colResultOther, force, hessian, meshId, vertexId, apply_friction);

		}
	}
}

void GAIA::VBDPhysics::updateCollisionInfoForVertex(TetMeshFEM* pMesh, IdType meshId, IdType vertexId)
{
	IdType surfaceVertexId = pMesh->tetVertIndicesToSurfaceVertIndices()(vertexId);

	if (surfaceVertexId > 0)
	{
		VBDCollisionDetectionResult& colResult = getCollisionDetectionResultFromSurfaceMeshId(meshId, surfaceVertexId);
		updateCollisionInfo(colResult);
	}

	CollisionRelationList& colRelations = getCollisionRelationList(meshId, vertexId);
	for (size_t iCollision = 0; iCollision < colRelations.size(); iCollision++)
	{
		CollisionRelation& colRelation = colRelations[iCollision];

		VBDCollisionDetectionResult& colResultOther = getCollisionDetectionResultFromSurfaceMeshId(colRelation.meshId, colRelation.surfaceVertexId);
		updateCollisionInfo(colResultOther);

	}
}

void GAIA::VBDPhysics::accumlateCollisionForceAndHessianPerCollision(int collisionId, int collisionVertexOrder, VBDCollisionDetectionResult& collisionResult,
	Vec3& force, Mat3& hessian, IdType meshId, IdType vertexId, bool apply_friction)
{
	CollidingPointInfo& collidingPt = collisionResult.collidingPts[collisionId];
	VBDCollisionInfo& collisionFH = collisionResult.collisionForceAndHessian[collisionId];

	if (!collidingPt.shortestPathFound)
	{
		return;
	}

	FloatingType b;
	CFloatingType k = physicsParams().collisionStiffness;
	CFloatingType collisionRadius = physicsParams().collisionAirDistance;
	const auto& n = collidingPt.closestPointNormal;
	const auto& diff = collisionFH.diff;
	const auto& frictionForce = collisionFH.frictionForce;
	const auto& frictionHessian = collisionFH.frictionHessian;

	FloatingType penetrationDepth = diff.dot(n);

	if (physicsParams().collisionEnergyType == 1)
	{
		penetrationDepth += physicsParams().collisionAirDistance;
	}

	if (penetrationDepth > 0)
	{
		switch (collisionVertexOrder)
		{
		case 0:
			b = -collidingPt.closestSurfacePtBarycentrics[0];
			break;
		case 1:
			b = -collidingPt.closestSurfacePtBarycentrics[1];
			break;
		case 2:
			b = -collidingPt.closestSurfacePtBarycentrics[2];
			break;
		case 3:
			b = 1.0f;
			break;
		default:
			break;
		}

		if (b != 0)
		{
			// Friction
			CFloatingType lambda = penetrationDepth * k;

			switch (physicsParams().collisionEnergyType)
			{
			case 0:
				force += k * b * diff;
				hessian += k * b * b * Mat3::Identity();
				break;
			case 1:
				// simplified point-plane
				// Penalty force
				debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
					collisionForce[meshId].col(vertexId) += k * b * penetrationDepth * n;
					});
				force += k * b * penetrationDepth * n;
				hessian += k * b * b * n * n.transpose();

				// Friction
				if (apply_friction) {
					debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
						this->frictionForce[meshId].col(vertexId) += frictionForce * (lambda * b);
						});
					force += frictionForce * (lambda * b);
					hessian += frictionHessian * (lambda * b * b);
				}
				break;
			default:
				break;
			}

		}
	}
}

VBDCollisionDetectionResult& GAIA::VBDPhysics::getCollisionDetectionResultFromSurfaceMeshId(int meshId, int surfaceVertexId)
{
	return collisionResultsAll[meshId][surfaceVertexId];
}

VBDCollisionDetectionResult& GAIA::VBDPhysics::getCollisionDetectionResultFromTetMeshId(int meshId, int vertexId)
{
	int surfaceVertexId = tMeshes[meshId]->tetVertIndicesToSurfaceVertIndices()(vertexId);
	assert(surfaceVertexId >= 0);
	return collisionResultsAll[meshId][surfaceVertexId];

}



CollisionRelationList& GAIA::VBDPhysics::getCollisionRelationList(int meshId, int vertexId)
{
	return activeColllisionList.vertexCollisionRelations[meshId][vertexId];
}

std::shared_ptr<ObjectParamsList> GAIA::VBDPhysics::createObjectParamsList()
{
	return std::make_shared<ObjectParamsListVBD>();
}

std::shared_ptr<BasePhysicsParams> GAIA::VBDPhysics::createPhysicsParams()
{
	return std::make_shared<VBDPhysicsParameters>();
}

std::shared_ptr<CollisionDetectionParamters> GAIA::VBDPhysics::createCollisionParams()
{
	return std::make_shared<CollisionDetectionParamters>();
}

std::shared_ptr<RunningTimeStatistics> GAIA::VBDPhysics::createRunningTimeStatistics()
{
	return std::make_shared<RunningTimeStatistics>();
}

ObjectParams::SharedPtr GAIA::ObjectParamsListVBD::createObjectParam(const std::string& materialName)
{
	ObjectParams::SharedPtr pObjParams;
	if (materialName == "NeoHookean")
	{
		pObjParams = std::make_shared<ObjectParametersVBDNeoHookean>();

	}
	else if (materialName == "MassSpring") {
		// pObjParams = std::make_shared<ObjectParametersEBDMassSpring>();

		std::cout << "Warning!!! Material name: MassSpring" << " not implemented yet! Skipping this object!\n";
		assert(false);
	}
	else
	{
		std::cout << "Warning!!! Material name: " << materialName << " not recognized! Skipping this material!\n";
		/*std::cout << "Using default material NeoHookean instead!\n";
		objectParams.push_back(std::make_shared<ObjectParamsPBDNeoHookean>());*/
	}
	return pObjParams;
}

//void VBDPhysics::testGPUCollisionHandlingCode()
//{
//	VBDPhysicsDataGPU* pPhysicsData = &vbdPhysicsDataGPU_forCPUDebug;
//	Vec3 forceCPU, forceGPU;
//	Mat3 hessianCPU, hessianGPU;
//
//	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
//	{
//		VBDBaseTetMeshGPU* pTetMeshGPU = vbdPhysicsDataGPU_forCPUDebug.tetMeshes[iMesh];
//
//		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
//
//		for (size_t iV = 0; iV < pTetMesh->numVertices(); iV++)
//		{
//			forceCPU.setZero();
//			forceGPU.setZero();
//			hessianCPU.setZero();
//			hessianGPU.setZero();
//
//			if (pTetMesh->activeCollisionMask[iV])
//			{
//				updateCollisionInfoForVertex(pTetMesh, iMesh, iV);
//				accumlateCollisionForceAndHessian(pTetMesh, iMesh, iV, forceCPU, hessianCPU);
//
//				printCPUCollisionDataForVertex(iMesh, iV);
//
//				CollisionDataGPU& collisionDataGPU = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(iV);
//				if (collisionDataGPU.activeColliding)
//				{
//					updateCollisionInfoGPU(pPhysicsData, iMesh, iV);
//				}
//
//				if (collisionDataGPU.numCollisionRelations)
//				{
//					for (int32_t iCollision = 0; iCollision < collisionDataGPU.numCollisionRelations; iCollision++)
//					{
//						const CollisionRelationGPU& collisionRelation = collisionDataGPU.collisionRelations[iCollision];
//						const VBDBaseTetMeshGPU* pTetMeshGPUOther = pPhysicsData->tetMeshes[collisionRelation.meshId];
//						const CollisionDataGPU& collisionDataOther = pTetMeshGPUOther->getCollisionData(collisionRelation.collisionPrimitiveId);
//
//						updateCollisionInfoGPU(pPhysicsData, collisionRelation.meshId, collisionRelation.collisionPrimitiveId);
//					}
//				}
//				accumulateCollisionForceAndHessian(&vbdPhysicsDataGPU_forCPUDebug, pTetMeshGPU,
//					iV, forceGPU.data(), hessianGPU.data());
//				printVBDPhysicsDataGPUForVertex(pPhysicsData, iMesh, iV);
//
//				// updateCollisionInfoGPU()
//
//				compareCPUandGPU(forceCPU, hessianCPU, forceGPU, hessianGPU);
//
//			}
//			//std::cout << "---------------------------------\n";
//
//			accumlateBoundaryForceAndHessian(pTetMesh, iMesh, iV, forceCPU, hessianCPU);
//			accumulateBoundaryForceAndHessianGPU(pPhysicsData, pTetMeshGPU, iV, forceGPU.data(), hessianGPU.data());
//			compareCPUandGPU(forceCPU, hessianCPU, forceGPU, hessianGPU);
//
//			VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
//			if (!pMesh->fixedMask[iV])
//			{
//				pMesh->accumlateInertiaForceAndHessian(iV, forceCPU, hessianCPU);
//				accumulateInertiaForceAndHessian(pTetMeshGPU, iV, physicsParams().dtSqrReciprocal, forceGPU.data(), hessianGPU.data());
//				compareCPUandGPU(forceCPU, hessianCPU, forceGPU, hessianGPU);
//
//				pMesh->accumlateMaterialForceAndHessian(iV, forceCPU, hessianCPU);
//				accumulateMaterialForceAndHessianForVertex_NeoHookean((VBDTetMeshNeoHookeanGPU*)pTetMeshGPU, iV, forceGPU.data(), hessianGPU.data(), pPhysicsData->dt);
//				//compareCPUandGPU(forceCPU, hessianCPU, forceGPU, hessianGPU);
//			}
//
//			// std::cout << "forceCPU: " << forceCPU.transpose() << "\n";
//			// std::cout << "hessianCPU:\n" << hessianCPU << "\n";
//			// 
//			// std::cout << "forceGPU: " << forceGPU.transpose() << "\n";
//			// std::cout << "hessianGPU:\n" << hessianGPU << "\n";
//
//
//
//			// std::cout << "---------------------------------\n";
//
//		}
//	}
//
//}

void GAIA::VBDPhysics::prepareCollisionDataCPUAndGPU()
{
	// clear the collision relations
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++)
	{
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		for (size_t iVert = 0; iVert < pTetMesh->numVertices(); iVert++)
		{
			CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(iVert);
			collisionDataCPUBuffer.numCollisionRelations = 0;
			collisionDataCPUBuffer.activeColliding = 0;
		}
	}
	activeColllisionList.clear();

	// record the activate collisions
	for (int iMesh = 0; iMesh < tMeshes.size(); iMesh++) {
		VBDBaseTetMesh* pTetMesh = tMeshes[iMesh].get();
		pTetMesh->activeCollisionMask.setZero();
		std::vector<VBDCollisionDetectionResult>& dcdResults = collisionResultsAll[iMesh];

		// todo: make this parallel using atomic operation
		// loop through all the v-f collisions
		for (int iSurfaceV = 0; iSurfaceV < pTetMesh->surfaceVIds().size(); ++iSurfaceV)
		{
			VBDCollisionDetectionResult& collisionResult = dcdResults[iSurfaceV];
			int32_t surfaceVIdTetMesh = pTetMesh->surfaceVIds()(iSurfaceV);
			CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(surfaceVIdTetMesh);

			if (collisionResult.numIntersections())
			{
				activeColllisionList.addToActiveCollisionList(collisionResult, 1);
				int iIntersection = 0;
				CollidingPointInfo& collidingPt = collisionResult.collidingPts[iIntersection];
				if (collidingPt.shortestPathFound)
				{
					//std::cout << "v side: collision mask of " << iMesh << "," << surfaceVIdTetMesh << " set to true\n";
					pTetMesh->activeCollisionMask(surfaceVIdTetMesh) = true;

					const int meshId_intersecting = collidingPt.intersectedMeshId;
					const int surfaceFaceId = collidingPt.closestSurfaceFaceId;

					collisionDataCPUBuffer.activeColliding = 1;

					collisionDataCPUBuffer.closestSurfaceFaceId = surfaceFaceId;
					collisionDataCPUBuffer.intersectedTMeshId = meshId_intersecting;

					VBDBaseTetMesh* pTetMesh_intersecting = tMeshes[meshId_intersecting].get();

					for (size_t iSurfaceFaceV = 0; iSurfaceFaceV < 3; iSurfaceFaceV++)
					{
						IdType faceVId = pTetMesh_intersecting->surfaceFacesTetMeshVIds()(iSurfaceFaceV, surfaceFaceId);

						collisionDataCPUBuffer.closestSurfaceFaceVIds[iSurfaceFaceV] = faceVId;

						int numColRelations = pTetMesh_intersecting->pTetMeshSharedBase->numCollisionRelation(faceVId);

						if (numColRelations < GPU_COLLISION_RELATION_PREALLOCATION_SIZE)
						{

							//std::cout << "f side: collision mask of " << meshId_intersecting << "," << faceVId << " set to true\n";
							CollisionDataGPU& collisionRelationOtherSide =
								pTetMesh_intersecting->pTetMeshSharedBase->getCollisionDataCPUBuffer(faceVId);
							pTetMesh_intersecting->activeCollisionMask(faceVId) = true;

							pTetMesh_intersecting->pTetMeshSharedBase->numCollisionRelation(faceVId)++;
							CollisionRelationGPU& collisionRelation =
								pTetMesh_intersecting->pTetMeshSharedBase->getCollisionRelationCPUBuffer(faceVId, numColRelations);
							collisionRelation.collisionPrimitiveId = surfaceVIdTetMesh;
							collisionRelation.collisionType = 0;
							collisionRelation.collisionVertexOrder = iSurfaceFaceV;
							collisionRelation.meshId = iMesh;

							CollisionRelationList& collisionRelationCPU = getCollisionRelationList(meshId_intersecting, faceVId);
							collisionRelationCPU.emplace_back();
							collisionRelationCPU.back().meshId = iMesh;
							collisionRelationCPU.back().surfaceVertexId = iSurfaceV;
							collisionRelationCPU.back().collisionId = iIntersection;
							collisionRelationCPU.back().collisionType = 0;
							collisionRelationCPU.back().collisionVertexOrder = iSurfaceFaceV;
						}
					}
				}
			}
		}// surface verts

		for (size_t iVert = 0; iVert < pTetMesh->numVertices(); iVert++)
		{
			CollisionDataGPU& collisionDataCPUBuffer = pTetMesh->pTetMeshSharedBase->getCollisionDataCPUBuffer(iVert);
			if (collisionDataCPUBuffer.activeColliding)
			{
				printf("Vert %d of mesh %d collide with mesh: %d, face: %d, %d, %d\n", iVert, iMesh, collisionDataCPUBuffer.intersectedTMeshId,
					collisionDataCPUBuffer.closestSurfaceFaceVIds[0],
					collisionDataCPUBuffer.closestSurfaceFaceVIds[1],
					collisionDataCPUBuffer.closestSurfaceFaceVIds[2]);
			}
			/*int numCollisionRelations = collisionDataCPUBuffer.numCollisionRelations;
			if (numCollisionRelations)
			{
				printf("Vertex %d has Number Collision Relations: %d\n", iVert, numCollisionRelations);
				for (size_t iRelation = 0; iRelation < numCollisionRelations; iRelation++)
				{
					CollisionRelationGPU& collisionRelation =
						pTetMesh->pTetMeshSharedBase->getCollisionRelationCPUBuffer(iVert, iRelation);
					printf("collisionPrimitiveId: %d, collisionVertexOrder: %d \n", collisionRelation.collisionPrimitiveId,
						collisionRelation.collisionVertexOrder);
				}
			}*/
		}

	}// meshes
}

void GAIA::VBDPhysics::printCPUCollisionDataForVertex(IdType meshId, IdType vertexId)
{
	printf("------Print collision info for vertex %d of mesh %d------\n", meshId, vertexId);
	VBDBaseTetMesh* pMesh = tMeshes[meshId].get();
	IdType surfaceVertexId = pMesh->tetVertIndicesToSurfaceVertIndices()(vertexId);

	if (surfaceVertexId > 0)
	{
		printf("v side:\n");
		VBDCollisionDetectionResult& colResult = getCollisionDetectionResultFromSurfaceMeshId(meshId, surfaceVertexId);
		if (colResult.numIntersections())
		{
			printCPUCollisionData(colResult);
		}
	}

	CollisionRelationList& colRelations = getCollisionRelationList(meshId, vertexId);
	for (size_t iCollision = 0; iCollision < colRelations.size(); iCollision++)
	{
		printf("f side (%d of %d):\n", iCollision, colRelations.size());

		CollisionRelation& colRelation = colRelations[iCollision];

		VBDCollisionDetectionResult& colResultOther = getCollisionDetectionResultFromSurfaceMeshId(colRelation.meshId, colRelation.surfaceVertexId);
		printCPUCollisionData(colResultOther);

	}
}

void GAIA::VBDPhysics::printCPUCollisionData(VBDCollisionDetectionResult& colResult, int iIntersection)
{
	const CollidingPointInfo& collidingPt = colResult.collidingPts[iIntersection];
	if (collidingPt.shortestPathFound)
	{
		VBDBaseTetMesh* pIntersectedTM = tMeshes[collidingPt.intersectedMeshId].get();
		int closestFaceId = collidingPt.closestSurfaceFaceId;
		printf("    collision %d-(%d, %d, %d)\n", colResult.idVQuery,
			pIntersectedTM->surfaceFacesTetMeshVIds()(0, closestFaceId),
			pIntersectedTM->surfaceFacesTetMeshVIds()(1, closestFaceId),
			pIntersectedTM->surfaceFacesTetMeshVIds()(2, closestFaceId));

		printf("    closestSurfacePt (%f, %f, %f)\n", collidingPt.closestSurfacePt[0],
			collidingPt.closestSurfacePt[1], collidingPt.closestSurfacePt[2]);
		printf("    closestPointNormal (%f, %f, %f)\n", collidingPt.closestPointNormal[0],
			collidingPt.closestPointNormal[1], collidingPt.closestPointNormal[2]);
		printf("    closestSurfacePtBarycentrics (%f, %f, %f)\n", collidingPt.closestSurfacePtBarycentrics[0],
			collidingPt.closestSurfacePtBarycentrics[1], collidingPt.closestSurfacePtBarycentrics[2]);
	}
}

void GAIA::VBDPhysics::compareCPUandGPU(Vec3& forceCPU, Mat3& hessianCPU, Vec3& forceGPU, Mat3& hessianGPU)
{
	CFloatingType eps = 1e-3f;
	CFloatingType epsAbs = 1e-2f;
	FloatingType maxForceNorm = std::max(forceCPU.norm(), forceGPU.norm());
	if (maxForceNorm)
	{
		if ((forceCPU - forceGPU).norm() / maxForceNorm > eps)
		{
			std::cout << "warning!!! cpu gpu force is different:\n"
				<< "forceCPU: " << forceCPU.transpose() << "\n"
				<< "forceGPU: " << forceGPU.transpose() << "\n";
		}
		assert((forceCPU - forceGPU).norm() / maxForceNorm < eps
			|| (forceCPU - forceGPU).norm() < epsAbs);
	}
	Mat3 diff = (hessianCPU - hessianGPU);
	FloatingType maHessianNorm = std::max(hessianCPU.norm(), hessianGPU.norm());
	if (maHessianNorm)
	{
		if (diff.norm() / maHessianNorm > eps)
		{
			std::cout << "warning!!! cpu gpu hessian is different:\n"
				<< "hessianCPU:\n" << hessianCPU << "\n"
				<< "hessianCPU:\n" << hessianCPU << "\n";
		}
		assert(diff.norm() / maHessianNorm < eps
			|| diff.norm() < epsAbs);

	}
}

void GAIA::VBDPhysics::outputPosVel()
{
	const auto debugFolder = getDebugFolder();
	std::string prefix = debugFolder + "/" + std::to_string(frameId) + "_" + std::to_string(substep) + "_" + std::to_string(iIter);
	std::string posFileName = prefix + ".pos";
	std::string velFileName = prefix + ".vel";
	std::ofstream posFile(posFileName, std::ios::binary);
	std::ofstream velFile(velFileName, std::ios::binary);
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDBaseTetMesh* pTM = tMeshes[iMesh].get();
		TVerticesMat velocity = (pTM->vertices() - pTM->positionsPrev()) / physicsParams().dt;
		velFile.write((char*)velocity.data(), sizeof(FloatingType) * velocity.size());
		posFile.write((char*)pTM->vertices().data(), sizeof(FloatingType) * pTM->vertices().size());
	}



}

void GAIA::VBDPhysics::outputForces()
{
	const auto debugFolder = getDebugFolder();
	std::string prefix = debugFolder + "/" + std::to_string(frameId) + "_" + std::to_string(substep) + "_" + std::to_string(iIter);
	std::string materialFileName = prefix + ".mf";
	std::ofstream materialFile(materialFileName, std::ios::binary);
	std::string inertiaFileName = prefix + ".if";
	std::ofstream inertiaFile(inertiaFileName, std::ios::binary);
	std::string boundaryFrictionFileName = prefix + ".bdff";
	std::ofstream boundaryFrictionFile(boundaryFrictionFileName, std::ios::binary);
	std::string boundaryCollisionFileName = prefix + ".bdcf";
	std::ofstream boundaryCollisionFile(boundaryCollisionFileName, std::ios::binary);
	std::string frictionFileName = prefix + ".ff";
	std::ofstream frictionFile(frictionFileName, std::ios::binary);
	std::string collisionFileName = prefix + ".cf";
	std::ofstream collisionFile(collisionFileName, std::ios::binary);

	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		const auto& mf = MaterialForce[iMesh];
		materialFile.write((char*)mf.data(), sizeof(FloatingType) * mf.size());
		const auto& iff = InertiaForce[iMesh];
		inertiaFile.write((char*)iff.data(), sizeof(FloatingType) * iff.size());
		const auto& bdff = boundaryFrictionForce[iMesh];
		boundaryFrictionFile.write((char*)bdff.data(), sizeof(FloatingType) * bdff.size());
		const auto& bdcf = boundaryCollisionForce[iMesh];
		boundaryCollisionFile.write((char*)bdcf.data(), sizeof(FloatingType) * bdcf.size());
		const auto& ff = frictionForce[iMesh];
		frictionFile.write((char*)ff.data(), sizeof(FloatingType) * ff.size());
		const auto& cf = collisionForce[iMesh];
		collisionFile.write((char*)cf.data(), sizeof(FloatingType) * cf.size());
	}
}

void GAIA::VBDPhysics::clearForces()
{
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		MaterialForce[iMesh].setZero();
		InertiaForce[iMesh].setZero();
		boundaryFrictionForce[iMesh].setZero();
		boundaryCollisionForce[iMesh].setZero();
		frictionForce[iMesh].setZero();
		collisionForce[iMesh].setZero();
	}
}

void GAIA::VBDPhysics::evaluateConvergence()
{
	NFloatingType e = 0.f;
	NFloatingType eInertia = 0.f;
	NFloatingType eElastic = 0.f;

	e = evaluateMeritEnergy(eInertia, eElastic);
	FloatingType avgForceNorm = 0.f;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iVert) {
			Vec3 force;
			Mat3 hessian;
			force = Vec3::Zero();
			if (!pMesh->fixedMask(iVert)) {
				pMesh->accumlateMaterialForceAndHessian2(iVert, force, hessian);

				force += pMesh->vertexMass(iVert) * (pMesh->inertia.col(iVert) - pMesh->vertex(iVert)) * (physicsParams().dtSqrReciprocal);
			}
			pMesh->vertexInternalForces.col(iVert) = force;
			});
		avgForceNorm += pMesh->vertexInternalForces.colwise().norm().sum();
		//for (size_t iVert = 0; iVert < pMesh->numVertices(); iVert++)
		//{
		//	avgForceNorm += pMesh->vertexInternalForces.col(iVert).norm();
		//}
	}

	avgForceNorm /= numAllVertices;

	std::cout << "Frame " << frameId << " step " << substep << " iteration: " << iIter
		<< ": energy " << e << " | inertia: " << eInertia << " | elastic: " << eElastic
		<< " | avg force norm: " << avgForceNorm
		<< std::endl;
	if (physicsParams().saveConvergenceEvaluationResults)
	{
		std::string outFolder = getDebugFolder();
		std::ostringstream aSs;
		aSs << "ConvergenceEvaluation_Frame" << std::setfill('0') << std::setw(8) << frameId << "_iter" << std::setw(4) << iIter;
		std::string outFile = outFolder + "/" + aSs.str() + ".json";

		ConvergenceStats stats(numTetMeshes());
		stats.energy[0] = e;
		stats.energy_elastic[0] = eElastic;
		stats.energy_inertia[0] = eInertia;
		stats.avgGradNorm[0] = avgForceNorm;
		stats.writeToJsonFile(outFile);
	}
}

void GAIA::VBDPhysics::evaluateConvergenceGPU()
{
	FloatingType e = 0.f;
	FloatingType eInertia = 0.f;
	FloatingType eElastic = 0.f;

	e = evaluateMeritEnergyGPU(eInertia, eElastic);
	FloatingType avgForceNorm = 0.f;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iVert) {
			Vec3 force;
			Mat3 hessian;
			force = Vec3::Zero();
			if (!pMesh->fixedMask(iVert)) {
				pMesh->accumlateMaterialForceAndHessian2(iVert, force, hessian);

				force += pMesh->vertexMass(iVert) * (pMesh->inertia.col(iVert) - pMesh->vertex(iVert)) * (physicsParams().dtSqrReciprocal);
			}
			pMesh->vertexInternalForces.col(iVert) = force;
			});
		avgForceNorm += pMesh->vertexInternalForces.colwise().norm().sum();
	}

	avgForceNorm /= numAllVertices;

	std::cout << "Frame " << frameId << " step " << substep << " iteration: " << iIter
		<< ": energy " << e << " | inertia: " << eInertia << " | elastic: " << eElastic
		<< " | avg force norm: " << avgForceNorm
		<< std::endl;

	if (physicsParams().saveConvergenceEvaluationResults)
	{
		std::string outFolder = getDebugFolder();
		std::ostringstream aSs;
		aSs << "ConvergenceEvaluation_Frame" << std::setfill('0') << std::setw(8) << frameId << "_iter" << std::setw(4) << iIter;
		std::string outFile = outFolder + "/" + aSs.str() + ".json";

		ConvergenceStats stats(numTetMeshes());
		stats.energy[0] = e;
		stats.energy_elastic[0] = eElastic;
		stats.energy_inertia[0] = eInertia;
		stats.avgGradNorm[0] = avgForceNorm;
		stats.writeToJsonFile(outFile);
	}
}

FloatingType GAIA::VBDPhysics::getAcceleratorOmega(int order, CFloatingType pho, CFloatingType prevOmega)
{
	switch (order)
	{
	case 1:
		return 1;
		break;
	case 2:
		return  2 / (2 - SQR(pho));
		break;
	default:
		assert(order > 0);

		return 4.f / (4 - SQR(pho) * prevOmega);;
		break;
	}
}

void GAIA::VBDPhysics::recordInitialPositionForAccelerator(bool sync)
{
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		// write the intial position to the previousPositionsForLineSearch for the first line search
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshShared->positionsPrevPrevIterBuffer->getGPUBuffer(), pMesh->pTetMeshSharedBase->vertPosBuffer->getGPUBuffer(),
			pMesh->pTetMeshShared->positionsPrevPrevIterBuffer->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
		// not all positionsNew will be modified, here we first fill it with the orginal position
		CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshShared->positionsNewBuffer->getGPUBuffer(), pMesh->pTetMeshSharedBase->vertPosBuffer->getGPUBuffer(),
			pMesh->pTetMeshShared->positionsPrevPrevIterBuffer->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
	}

	if (sync)
	{
		CHECK_CUDA_ERROR(cudaDeviceSynchronize());
	}
}

void GAIA::VBDPhysics::recordPrevIterPositionsAccelerator(bool sync)
{
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshShared->positionsPrevIterBuffer->getGPUBuffer(), pMesh->pTetMeshShared->vertPosBuffer->getGPUBuffer(),
			pMesh->pTetMeshShared->positionsPrevIterBuffer->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
	}

	if (sync)
	{
		CHECK_CUDA_ERROR(cudaDeviceSynchronize());
	}
}

void GAIA::VBDPhysics::computeElasticForceHessian()
{
	auto energyIter = &pNewtonAssembler->elasticEnergy[0];
	auto forceIter = &pNewtonAssembler->elasticForce[0];
	auto hessianIter = &pNewtonAssembler->elasticHessian[0];
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		//cpu_parallel_for(0, pMesh->numTets(), [&](int iTet) {
		//	pMesh->validateElasticGradientHessianF(iTet);
		//	});
		//cpu_parallel_for(0, pMesh->numTets(), [&](int iTet) {
		//	pMesh->validateElasticForceHessian(iTet);
		//	});
#ifdef USE_DOUBLE
		cpu_parallel_for(0, pMesh->numTets(), [&](int iTet) {
			pMesh->computeElasticForceHessianDouble(iTet, *(energyIter + iTet), *(forceIter + iTet), *(hessianIter + iTet));
			});
#else
		cpu_parallel_for(0, pMesh->numTets(), [&](int iTet) {
			pMesh->computeElasticForceHessian(iTet, *(energyIter + iTet), *(forceIter + iTet), *(hessianIter + iTet));
			});
#endif
		//for (size_t iTet = 0; iTet < pMesh->numTets(); iTet++)
		//{
		//	pMesh->computeElasticForceHessian(iTet, *forceIter, *hessianIter);
		//	forceIter++;
		//	hessianIter++;
		//}

		energyIter += pMesh->numTets();
		forceIter += pMesh->numTets();
		hessianIter += pMesh->numTets();
	}
}

void GAIA::VBDPhysics::computeElasticEnergy()
{
	auto energyIter = &pNewtonAssembler->elasticEnergy[0];
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numTets(), [&](int iTet) {
			pMesh->computeElasticEnergy(iTet, *(energyIter + iTet));
			});
		energyIter += pMesh->numTets();
	}
}

void GAIA::VBDPhysics::fillNewtonSystem()
{
	pNewtonAssembler->newtonForce.setZero();
	fillNewtonForce();

	auto data_ptr = pNewtonAssembler->newtonHessianAll.valuePtr();
	memset(data_ptr, 0, pNewtonAssembler->newtonHessianAll.nonZeros() * sizeof(FloatingType));
	fillNewtonHessianDiagonal();
	fillNewtonHessianOffDiagonal();
	// std::cout << newtonHessian.block(0, 0, 9, 9) << "\n";
}

void GAIA::VBDPhysics::fillNewtonForce()
{
	NVec12* elasticForcePtr = &pNewtonAssembler->elasticForce[0];
	NFloatingType* forcePtr = pNewtonAssembler->newtonForce.data();
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {
			if (!pMesh->fixedMask[iV])
			{
				Eigen::Map<NVec3> force(forcePtr + iV * 3);
				// External force and inertia force
				force = (pMesh->vertexExternalForces.col(iV) + pMesh->vertexMass(iV) * (pMesh->inertia.col(iV) - pMesh->vertex(iV)) * (physicsParams().dtSqrReciprocal)).cast<NFloatingType>();
				// Elastic force
				const size_t numNeiTest = pMesh->getNumVertexNeighborTets(iV);
				for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++) {
					const auto tetId = pMesh->getVertexNeighborTet(iV, iNeiTet);
					const auto vertedTetVId = pMesh->getVertexNeighborTetVertexOrder(iV, iNeiTet);
					force += elasticForcePtr[tetId].segment<3>(vertedTetVId * 3);
				}
			}
			});
		elasticForcePtr += pMesh->numTets();
		forcePtr += pMesh->numVertices() * 3;
	}
}

void GAIA::VBDPhysics::fillNewtonHessianDiagonal()
{
	NMat12* elasticHessianPtr = &pNewtonAssembler->elasticHessian[0];
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numVertices(), [&](int iV) {

			NMat3 hessian = NMat3::Identity() * pMesh->vertexMass(iV) * (physicsParams().dtSqrReciprocal);
			if (!pMesh->fixedMask[iV])
			{
				// Elastic hessian
				const size_t numNeiTets = pMesh->getNumVertexNeighborTets(iV);
				for (size_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++) {
					const auto tetId = pMesh->getVertexNeighborTet(iV, iNeiTet);
					const auto vertedTetVId = pMesh->getVertexNeighborTetVertexOrder(iV, iNeiTet);
					hessian += elasticHessianPtr[tetId].block<3, 3>(vertedTetVId * 3, vertedTetVId * 3);
				}
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
		elasticHessianPtr += pMesh->numTets();
	}
}

void GAIA::VBDPhysics::fillNewtonHessianOffDiagonal()
{
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		cpu_parallel_for(0, pMesh->numEdges(), [&](int iE) {
			int iV = pMesh->edges()(0, iE);
			int jV = pMesh->edges()(1, iE);
			if (!pMesh->fixedMask[iV] && !pMesh->fixedMask[jV])
			{
				NMat3 hessian = NMat3::Zero();
				// Elastic hessian
				const size_t numNeiTets = pMesh->getNumEdgeNeighborTets(iE);
				for (size_t iNeiTet = 0; iNeiTet < numNeiTets; iNeiTet++) {
					const auto tetId = pMesh->getEdgeNeighborTet(iE, iNeiTet);
					int corner0, corner1;
					pMesh->getEdgeNeighborTetVertexOrder(iE, iNeiTet, corner0, corner1);
					hessian += pNewtonAssembler->elasticHessian[tetId].block<3, 3>(corner0 * 3, corner1 * 3);
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

void GAIA::VBDPhysics::updatePositions(const VecDynamic& dx)
{
	int offset = 0;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		pMesh->vertices() += physicsParams().stepSize 
			* Eigen::Map<const TVerticesMat>(dx.data() + offset, 3, pMesh->numVertices());
		offset += pMesh->numVertices() * 3;
	}
}

GAIA::NFloatingType GAIA::VBDPhysics::evaluateMeritEnergy(NFloatingType& eInertia, NFloatingType& eElastic, bool elasticReady)
{
	if (!elasticReady)computeElasticEnergy();
	eElastic = 0;
	for (auto& e : pNewtonAssembler->elasticEnergy)
	{
		eElastic += e;
	}
	eInertia = 0;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		eInertia += pMesh->computeInertiaEnergy();
	}
	return eElastic + eInertia;
}

GAIA::NFloatingType GAIA::VBDPhysics::newtonLineSearch(const VecDynamic& dx, NFloatingType E0, FloatingType alpha,
	FloatingType c, FloatingType tau, int maxNumIters, FloatingType& stepSizeOut)
{
	FloatingType m = dx.squaredNorm();

	std::vector<TVerticesMat> orgPos{};
	orgPos.reserve(numTetMeshes());
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		orgPos.push_back(tMeshes[iMesh]->vertices());
	}

	NFloatingType e = E0;
	NFloatingType eInertia = 0;
	NFloatingType eElastic = 0;

	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
		std::cout << "Initial Energy: " << E0 << std::endl;
		});
	for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
	{
		int offset = 0;
		for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
		{
			VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
			pMesh->vertices() = orgPos[iMesh] + alpha * Eigen::Map<const TVerticesMat>(dx.data() + offset, 3, pMesh->numVertices());
			offset += pMesh->numVertices() * 3;
		}
		e = evaluateMeritEnergy(eInertia, eElastic);
		debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
			std::cout << "alpha: " << alpha << ", energy: " << e << ", inertia: " << eInertia << ", elastic: " << eElastic << std::endl;
			});

		// first Wolfe condition 
		if (e < E0 - alpha * c * m)
		{
			// std::cout << "step size for vertex " << vId << ": " << alpha << "\n";
			break;
		}
		else
		{
			alpha = alpha * tau;
		}
	}

	stepSizeOut = alpha;

	return e;
}

FloatingType GAIA::VBDPhysics::evaluateMeritEnergyGPU(FloatingType& eInertia, FloatingType& eElastic, bool elasticReady)
{
	evaluateElasticEnergyGPU(pPhysicsDataGPUBuffer->getData(), pLineSearchUtilities->tetAllParallelGroupsBuffer->getGPUBuffer(), numAllTets, pLineSearchUtilities->tetElasticEnergyBuffer->getGPUBuffer(),
		physicsParams().numThreadsVBDSolve, cudaStream);
	pLineSearchUtilities->tetElasticEnergyBuffer->toCPU(false, cudaStream);
	syncAllToCPUVertPosOnly(true);

	eInertia = 0.f;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		eInertia += pMesh->computeInertiaEnergy();
	}

	eElastic = 0.f;
	for (size_t iTet = 0; iTet < numAllTets; iTet++)
	{
		eElastic += pLineSearchUtilities->tetElasticEnergyBuffer->getCPUBuffer()[iTet];
	}
	return eInertia + eElastic;
}

void GAIA::VBDPhysics::GDBackupPositions(bool sync, CFloatingType omega)
{
	pGDSolverUtilities->omegaPrev = omega;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		// write the intial position to the previousPositionsForLineSearch for the first line search
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		CHECK_CUDA_ERROR(cudaMemcpyAsync(pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->getGPUBuffer(), pMesh->pTetMeshSharedBase->vertPosBuffer->getGPUBuffer(),
			pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));

		if (physicsParams().useAccelerator)
		{
			CHECK_CUDA_ERROR(cudaMemcpyAsync(pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->getGPUBuffer(), pMesh->pTetMeshSharedBase->positionsPrevPrevIterBuffer->getGPUBuffer(),
				pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
		}
	}

	if (sync)
	{
		CHECK_CUDA_ERROR(cudaDeviceSynchronize());
	}

	//syncAllToCPUVertPosOnly(true);
	//for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	//{
	//	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
	//	//pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->toCPU(true, cudaStream);
	//	//Eigen::Map<Eigen::Matrix<FloatingType, 3, -1>> prevpos(pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->getCPUBuffer(), 3, pMesh->numVertices());

	//	//std::cout << prevpos.transpose();
	//	//std::cout << (pMesh->vertices() - prevpos).transpose();

	//	pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->toCPU(true, cudaStream);
	//	Eigen::Map<Eigen::Matrix<FloatingType, 3, -1>> posPrevIter(pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->getCPUBuffer(), 3, pMesh->numVertices());
	//	std::cout << posPrevIter.transpose();
	//	std::cout << (pMesh->vertices() - posPrevIter).transpose();
	//}



}

void GAIA::VBDPhysics::GDRevertToBackupPositions(bool sync, FloatingType& omega)
{
	omega = pGDSolverUtilities->omegaPrev;
	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	{
		// write the intial position to the previousPositionsForLineSearch for the first line search
		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
		CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshSharedBase->vertPosBuffer->getGPUBuffer(), pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->getGPUBuffer(),
			pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
		if (physicsParams().useAccelerator)
		{
			CHECK_CUDA_ERROR(cudaMemcpyAsync(pMesh->pTetMeshSharedBase->positionsPrevPrevIterBuffer->getGPUBuffer(), pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->getGPUBuffer(),
				pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->nBytes(), cudaMemcpyDeviceToDevice, cudaStream));
		}
	}
	if (sync)
	{
		CHECK_CUDA_ERROR(cudaDeviceSynchronize());
	}

	//syncAllToCPUVertPosOnly(true);
	//for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	//{
	//	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
	//	//pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->toCPU(true, cudaStream);
	//	//Eigen::Map<Eigen::Matrix<FloatingType, 3, -1>> prevpos(pGDSolverUtilities->previousPositionsForLineSearchBuffer[iMesh]->getCPUBuffer(), 3, pMesh->numVertices());

	//	//std::cout << prevpos.transpose();
	//	//std::cout << (pMesh->vertices() - prevpos).transpose();

	//	pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->toCPU(true, cudaStream);
	//	pMesh->pTetMeshShared->positionsPrevPrevIterBuffer->toCPU(true);
	//	Eigen::Map<Eigen::Matrix<FloatingType, 3, -1>> posPrevIterInMesh(pMesh->pTetMeshShared->positionsPrevPrevIterBuffer->getCPUBuffer(), 3, pMesh->numVertices());


	//	Eigen::Map<Eigen::Matrix<FloatingType, 3, -1>> posPrevIter(pGDSolverUtilities->prevPrevPositionsForAcceleratorBuffer[iMesh]->getCPUBuffer(), 3, pMesh->numVertices());
	//	std::cout << posPrevIter.transpose();
	//	std::cout << (posPrevIterInMesh - posPrevIter).transpose();
	//}

}

FloatingType GAIA::VBDPhysics::GDLineSearch(GD_SolverUtilities::SharedPtr pGDSolverUtilities, FloatingType& E0, FloatingType& EInertia_prev, FloatingType& EElastic_prev,
	FloatingType alpha, FloatingType c, FloatingType tau, int maxNumIters)
{
	//FloatingType m = 0;
	//
	//for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	//{
	//	m += pGDSolverUtilities->dxs[iMesh].squaredNorm();
	//}

	//FloatingType e = E0;
	//FloatingType eInertia = 0;
	//FloatingType eElastic = 0;

	//debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
	//	std::cout << "Initial Energy: " << E0 << std::endl;
	//	});
	//for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
	//{
	//	for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	//	{
	//		VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
	//		pMesh->vertices() = pGDSolverUtilities->previousPositionsForLineSearch[iMesh] 
	//			+ alpha * pGDSolverUtilities->dxs[iMesh];
	//	}
	//	e = evaluateMeritEnergy(eInertia, eElastic);
	//	debugOperation(DEBUG_LVL_DEBUG_VEBOSE, [&]() {
	//		std::cout << "alpha: " << alpha << ", energy: " << e << ", inertia: " << eInertia << ", elastic: " << eElastic << std::endl;
	//		});

	//	// first Wolfe condition 
	//	if (e < E0 - alpha * c * m)
	//	{
	//		// std::cout << "step size for vertex " << vId << ": " << alpha << "\n";
	//		break;
	//	}
	//	else
	//	{
	//		alpha = alpha * tau;
	//	}
	//}

	//// write the solution to the previousPositionsForLineSearch for the next line search
	//for (size_t iMesh = 0; iMesh < numTetMeshes(); iMesh++)
	//{
	//	VBDTetMeshNeoHookean* pMesh = (VBDTetMeshNeoHookean*)tMeshes[iMesh].get();
	//	pGDSolverUtilities->previousPositionsForLineSearch[iMesh] = pMesh->vertices();
	//}

	//E0 = e;
	//EInertia_prev = eInertia;
	//EElastic_prev = eElastic;

	//return alpha;

	return 0.f;
}

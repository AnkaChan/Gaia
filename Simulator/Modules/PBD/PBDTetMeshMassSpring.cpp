#include "PBDTetMeshMassSpring.h"
#include "PBDTetMeshMassSpringCompute.h"

#include <random>
#include <chrono>

bool GAIA::ObjectParamsPBDMassSpring::fromJson(nlohmann::json& objectParam)
{
	ObjectParamsPBD::fromJson(objectParam);

	EXTRACT_FROM_JSON(objectParam, springCompliance);
	EXTRACT_FROM_JSON(objectParam, springDamping);

	return true;
}

bool GAIA::ObjectParamsPBDMassSpring::toJson(nlohmann::json& objectParam)
{
	ObjectParamsPBD::toJson(objectParam);

	PUT_TO_JSON(objectParam, springCompliance);
	PUT_TO_JSON(objectParam, springDamping);

	return true;
}

void GAIA::PBDTetMeshMassSpring::initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF, PBDPhysics* inPPBDPhysics)
{
	PBDTetMeshFEM::initialize(inMaterialParams, pTM_MF, inPPBDPhysics);
	pObjectParamsMaterial = std::static_pointer_cast<ObjectParamsPBDMassSpring>(pObjectParams);

	pObjectParams = inMaterialParams;

	springLambdas.resize(numEdges());
	orgLengths.resize(numEdges());

	for (size_t iEdge = 0; iEdge < numEdges(); iEdge++)
	{
		springLambdas(iEdge) = 0.f;
		orgLengths(iEdge) = (mVertPos.col(edges()(0, iEdge)) - mVertPos.col(edges()(1, iEdge))).norm();
	}

	// std::cout << edges.transpose() << std::endl;

	initializeGPUTetMesh(this, &tetMeshGPU);
	pTetMeshGPUBuffer = std::make_shared<DeviceClassBuffer<PBDTetMeshFEMGPUMassSpring>>();
	pTetMeshGPUBuffer->fromCPU(&tetMeshGPU);

	initializeGPUTetMeshOnCPU(this, &tetMeshGPU_forCPU);
}

void GAIA::PBDTetMeshMassSpring::applyInitialGuess()
{
	PBDTetMeshFEM::applyInitialGuess();
	initializeLambda();
}

void GAIA::PBDTetMeshMassSpring::syncToGPU(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToGPU(false, stream);
	springLambdasBuffer->toGPU(false, stream);
	if (sync) {
		cudaStreamSynchronize(stream);
	}

}

void GAIA::PBDTetMeshMassSpring::syncToCPU(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPU(sync, stream);
}

void GAIA::PBDTetMeshMassSpring::syncToCPUVertPosOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPUVertPosOnly(sync, stream);
}

void GAIA::PBDTetMeshMassSpring::syncToGPUVertPosOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToGPUVertPosOnly(sync, stream);
}

void GAIA::PBDTetMeshMassSpring::syncToCPUInversionSignOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPUInversionSignOnly(sync, stream);
}

void GAIA::PBDTetMeshMassSpring::solveMaterialGPUAggregated(int numThreads, cudaStream_t stream)
{
	solveSpringConstraintsAggregatedGPU(numThreads, stream, edgesColoringCategories().front().size(), getTetMeshGPU());
	solveVolumeConstraintsAggregatedGPU(numThreads, stream, tetsColoringCategories().front().size(), getTetMeshGPU());
}

void GAIA::PBDTetMeshMassSpring::solveMaterialConstraintGPU_ParallelizationGroup(int iGroup, int numThreads, cudaStream_t stream)
{
	if (iGroup < edgesColoringCategories().size())
		// solve edges first
	{
		solveSpringConstraintsGPU(edgesColoringEachCategoryBuffer[iGroup]->getGPUBuffer(), edgesColoringEachCategoryBuffer[iGroup]->getSize(),
			numThreads, getTetMeshGPU(), stream);
	}
	else 
		// solve volume constraints the second
	{
		iGroup = iGroup - edgesColoringCategories().size();
		solveVolumeConstraintsGPU(tetsColoringEachCategoryBuffer[iGroup]->getGPUBuffer(), tetsColoringEachCategoryBuffer[iGroup]->getSize(),
			numThreads, getTetMeshGPU(), stream);
	}

}

void GAIA::PBDTetMeshMassSpring::solveMaterialConstraintGPU_debugOnCPU(PBDTetMeshFEMGPU* pExternGPUTMeshOnCPU)
{
	for (size_t iParalGroup = 0; iParalGroup < edgesColoringCategories().size(); iParalGroup++)
	{
		solveSpringConstraintsCPU(edgesColoringCategories()[iParalGroup].data(), edgesColoringCategories()[iParalGroup].size(), tetMeshGPU_forCPU);
	}

	for (size_t iParalGroup = 0; iParalGroup < tetsColoringCategories().size(); iParalGroup++)
	{
		solveVolumeConstraintsCPU(tetsColoringCategories()[iParalGroup].data(), tetsColoringCategories()[iParalGroup].size(), tetMeshGPU_forCPU);
	}
}

size_t GAIA::PBDTetMeshMassSpring::numAllParallelizationGroups()
{
	return edgesColoringEachCategoryBuffer.size() + tetsColoringEachCategoryBuffer.size();
}

size_t GAIA::PBDTetMeshMassSpring::numTetParallelizationGroups()
{
	return PBDTetMeshFEMShared::numTetParallelizationGroups();
}

size_t GAIA::PBDTetMeshMassSpring::getTetParallelizationGroupSize(int32_t iGroup)
{
	return PBDTetMeshFEMShared::getTetParallelizationGroupSize(iGroup);
}

int32_t* GAIA::PBDTetMeshMassSpring::getTetParallelizationGroupBufferGPU(int32_t iGroup)
{
	return  PBDTetMeshFEMShared::getTetParallelizationGroupBufferGPU(iGroup);
}

GAIA::FloatingType GAIA::PBDTetMeshMassSpring::evaluateElasticEnergy()
{
	FloatingType energy = 0.f;
	for (size_t iEdge = 0; iEdge < edges().cols(); iEdge++)
	{
		IdType v1 = edges()(0, iEdge);
		IdType v2 = edges()(1, iEdge);
		Vec3 diff = vertex(v1) - vertex(v2);
		FloatingType l = diff.norm();
		FloatingType l0 = orgLengths(iEdge);

		energy += 0.5 * (1.f / pObjectParamsMaterial->springCompliance) * (l - l0) * (l - l0);
	}

	//for (size_t iTet = 0; iTet < numTets(); iTet++)
	//{
	//	int id0 = mTetVIds(0, iTet);
	//	int id1 = mTetVIds(1, iTet);
	//	int id2 = mTetVIds(2, iTet);
	//	int id3 = mTetVIds(3, iTet);

	//	Vec3 diff1 = vertex(id1) - vertex(id0);
	//	Vec3 diff2 = vertex(id2) - vertex(id0);
	//	Vec3 diff3 = vertex(id3) - vertex(id0);

	//	Eigen::Matrix<FloatingType, 3, 3> Dm;
	//	Dm << diff1[0], diff2[0], diff3[0],
	//		diff1[1], diff2[1], diff3[1],
	//		diff1[2], diff2[2], diff3[2]
	//		;

	//	int indToDsInv = iTet * 9;
	//	Eigen::Matrix<FloatingType, 3, 3> DsInv(DmInvs.data() + indToDsInv);
	//	Eigen::Matrix<FloatingType, 3, 3> F = Dm * DsInv;

	//	// Neo-Hookean Models
	//	FloatingType PHI_H = 0.f;
	//	if (pObjectParamsMaterial->volCompliance != 0)
	//	{
	//		PHI_H = (1 / pObjectParamsMaterial->volCompliance) * pow((F.determinant() - 1), 2);
	//	}
	//	energy += (PHI_H) * Dm.determinant() / 6;
	//}

	return energy;
}

void GAIA::PBDTetMeshMassSpring::initializeLambda()
{
	springLambdas = decltype(springLambdas)::Zero(springLambdas.size());
}

GAIA::ObjectParamsPBDMassSpring::SharedPtr GAIA::PBDTetMeshMassSpring::getObjectParams()
{
	return pObjectParamsMaterial;
}

void GAIA::PBDTetMeshFEMShared_MassSpring::copyAllToGPU()
{
	PBDTetMeshFEMShared::copyAllToGPU();
	for (size_t iCategory = 0; iCategory < edgesColoringEachCategoryBuffer.size(); iCategory++)
	{
		edgesColoringEachCategoryBuffer[iCategory]->toGPU(false);

	}
	springLambdasBuffer->toGPU(false);
	orgLengthsBuffer->toGPU(false);
	edgesBuffer->toGPU(false);
	edgesColoringCategoriesPointersBuffer->toGPU(false);
	edgesColoringEachCategorySizeBuffer->toGPU(true);
}

void GAIA::PBDTetMeshFEMShared_MassSpring::initializeGPUTetMesh(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	PBDTetMeshFEMShared::initializeGPUTetMesh((PBDTetMeshFEM*)pTetMesh, pTetMeshGPU);
	std::cout << "Initializing GPU tet mesh for MassSpring material.\n";

	edgesBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->edges().size(), true, pTetMesh->edges().data());
	edgesBuffer->toGPU();

	springLambdasBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->springLambdas.size(), true, pTetMesh->springLambdas.data());
	springLambdasBuffer->toGPU();

	orgLengthsBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->orgLengths.size(), true, pTetMesh->orgLengths.data());
	orgLengthsBuffer->toGPU();

	edgesColoringCategoriesPointersBuffer = std::make_shared<ManagedBuffer<int32_t*>>(pTetMesh->edgesColoringCategories().size(), true);
	edgesColoringEachCategorySizeBuffer = std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->edgesColoringCategories().size(), true);

	for (size_t iColorCategory = 0; iColorCategory < pTetMesh->edgesColoringCategories().size(); iColorCategory++)
	{
		edgesColoringEachCategoryBuffer.push_back(
			std::make_shared<ManagedBuffer<int32_t>>(pTetMesh->edgesColoringCategories()[iColorCategory].size(), true, pTetMesh->edgesColoringCategories()[iColorCategory].data())
		);
		// copy CPU coloring buffer to GPU
		edgesColoringEachCategoryBuffer.back()->toGPU();

		edgesColoringCategoriesPointersBuffer->getCPUBuffer()[iColorCategory] = edgesColoringEachCategoryBuffer.back()->getGPUBuffer();
		edgesColoringEachCategorySizeBuffer->getCPUBuffer()[iColorCategory] = edgesColoringEachCategoryBuffer.back()->getSize();
	}
	edgesColoringCategoriesPointersBuffer->toGPU();
	edgesColoringEachCategorySizeBuffer->toGPU();

	SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, edges);
	SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, orgLengths);
	SET_GPU_MESH_BUFER_PTR(pTetMeshGPU, springLambdas);

	pTetMeshGPU->nEdges = pTetMesh->numEdges();
	pTetMeshGPU->springCompliance = pTetMesh->getObjectParams()->springCompliance;
	pTetMeshGPU->springDamping = pTetMesh->getObjectParams()->springDamping;
	pTetMeshGPU->dt = pTetMesh->dt;

	pTetMeshGPU->numEdgesColoredCatergories = pTetMesh->edgesColoringCategories().size();
	pTetMeshGPU->edgesColoringCategoriesPointers = edgesColoringCategoriesPointersBuffer->getGPUBuffer();
	pTetMeshGPU->edgesColoringEachCategorySize = edgesColoringEachCategorySizeBuffer->getGPUBuffer();
}

void GAIA::PBDTetMeshFEMShared_MassSpring::initializeGPUTetMeshOnCPU(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
{
	PBDTetMeshFEMShared::initializeGPUTetMeshForDebug((PBDTetMeshFEM*)pTetMesh, pTetMeshGPU);

	pTetMeshGPU->edges = pTetMesh->edges().data();
	pTetMeshGPU->springLambdas = pTetMesh->springLambdas.data();
	pTetMeshGPU->orgLengths = pTetMesh->orgLengths.data();

	pTetMeshGPU->nEdges = pTetMesh->numEdges();
	pTetMeshGPU->springCompliance = pTetMesh->getObjectParams()->springCompliance;
	pTetMeshGPU->springDamping = pTetMesh->getObjectParams()->springDamping;
	pTetMeshGPU->dt = pTetMesh->dt;
}

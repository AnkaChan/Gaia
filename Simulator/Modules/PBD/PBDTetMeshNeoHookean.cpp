#pragma once

#include "PBDTetMeshNeoHookean.h"
#include "PBDTetMeshNeoHookeanCompute.h"
#include "PBDTetMeshGeneralCompute.h"
#include <CuMatrix/MatrixOps/CuMatrix.h>
#include <CuMatrix/Buffers/ManagedBuffer.h>
#include <memory>
#include "MeshFrame/Parser/Parser.h"


using namespace GAIA;

bool GAIA::ObjectParamsPBDNeoHookean::fromJson(nlohmann::json& objectParam)
{
	ObjectParamsPBD::fromJson(objectParam);

	EXTRACT_FROM_JSON(objectParam, devCompliance);
	EXTRACT_FROM_JSON(objectParam, volCompliance);
	EXTRACT_FROM_JSON(objectParam, devDamping);
	EXTRACT_FROM_JSON(objectParam, volDamping);
	EXTRACT_FROM_JSON(objectParam, tetsColoringCategoriesPath);

	return true;
}

bool GAIA::ObjectParamsPBDNeoHookean::toJson(nlohmann::json& objectParam)
{
	ObjectParamsPBD::toJson(objectParam);

	PUT_TO_JSON(objectParam, devCompliance);
	PUT_TO_JSON(objectParam, volCompliance);
	PUT_TO_JSON(objectParam, devDamping);
	PUT_TO_JSON(objectParam, volDamping);
	PUT_TO_JSON(objectParam, tetsColoringCategoriesPath);

	return true;
}

void GAIA::PBDTetMeshNeoHookean::initialize(ObjectParams::SharedPtr inMaterialParams, std::shared_ptr<TetMeshMF> pTM_MF, PBDPhysics* inPPBDPhysics)
{
	PBDTetMeshFEM::initialize(inMaterialParams, pTM_MF, inPPBDPhysics);
	pObjectParams = inMaterialParams;
	pObjectParamsNeoHookean = std::static_pointer_cast<ObjectParamsPBDNeoHookean>(pObjectParams);

	devLambdas.resize(m_nTets);
	volLambdas.resize(m_nTets);

	initializeGPUTetMesh(this, &tetMeshGPU);
	pTetMeshGPUBuffer = std::make_shared<DeviceClassBuffer<TetMeshFEMGPU_NeoHookean>>();
	pTetMeshGPUBuffer->fromCPU(&tetMeshGPU);

	initializeGPUTetMeshForDebug(this, &tetMeshGPU_forCPU);
}

void GAIA::PBDTetMeshNeoHookean::applyInitialGuess()
{
	PBDTetMeshFEM::applyInitialGuess();

	initializeLambda();
}

void GAIA::PBDTetMeshNeoHookean::syncToGPU(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToGPU(false, stream);
	devLambdasBuffer->toGPU(false, stream);
	volLambdasBuffer->toGPU(false, stream);

	if (sync)
	{
		cudaStreamSynchronize(stream);
	}
}

void GAIA::PBDTetMeshNeoHookean::syncToCPU(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPU(sync, stream);
#ifdef EVAL_TET_ACCESS_COUNT
	tetAccessCountBuffer->toCPU();
	for (int iTet = 0; iTet < tetAccessCountBuffer->getSize() - 1; ++iTet) {
		//printf("%d\n", tetAccessCountBuffer->getCPUBuffer()[iTet]);
		assert(tetAccessCountBuffer->getCPUBuffer()[iTet] == tetAccessCountBuffer->getCPUBuffer()[iTet + 1]);
	}
	vertAccessCountBuffer->toCPU();
	for (int iVert = 0; iVert < vertAccessCountBuffer->getSize() - 1; ++iVert) {
		int numNeiT = 0;
		numNeiT+=this->idVertex(iVert)->tvertices()->size();
		//printf("access count: %d : num NeiT: %d\n", vertAccessCountBuffer->getCPUBuffer()[iVert], numNeiT);
		assert(vertAccessCountBuffer->getCPUBuffer()[iVert] == numNeiT * tetAccessCountBuffer->getCPUBuffer()[0]);
	}
#endif // !EVAL_TET_ACCESS_COUNT
}

void GAIA::PBDTetMeshNeoHookean::syncToCPUVertPosOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPUVertPosOnly(sync, stream);
}

void GAIA::PBDTetMeshNeoHookean::syncToGPUVertPosOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToGPUVertPosOnly(sync, stream);

}

void GAIA::PBDTetMeshNeoHookean::syncToCPUInversionSignOnly(bool sync, cudaStream_t stream)
{
	PBDTetMeshFEMShared::syncToCPUInversionSignOnly(sync, stream);
}


void GAIA::PBDTetMeshNeoHookean::solveMaterialGPUAggregated(int numThreads, cudaStream_t stream)
{
	solveNeoHookeanMaterialAggregatedGPU(numThreads, stream, tetsColoringCategories().front().size(), getTetMeshGPU());
}

void GAIA::PBDTetMeshNeoHookean::solveMaterialConstraintGPU_ParallelizationGroup(int iGroup, int numThreads, cudaStream_t stream)
{	
	// v1
	solveMaterial(tetsColoringEachCategoryBuffer[iGroup]->getGPUBuffer(), tetsColoringEachCategoryBuffer[iGroup]->getSize(),
		numThreads, getTetMeshGPU(), stream);
	//solveMaterial(tetsColoringEachCategoryBuffer[iGroup]->getGPUBuffer(), tetsColoringEachCategoryBuffer[iGroup]->getSize(), 
	//		numThreads, &tetMeshGPU, stream);


	// for debugging:

	//if (iGroup == 0) {
	//	printf("Size of TetMeshFEMGPU_NeoHookean on CPU (.cpp): %d\n", sizeof(TetMeshFEMGPU_NeoHookean));
	//	solveMaterial_oneThreadForLoop(tetsColoringEachCategoryBuffer[iGroup]->getGPUBuffer(), tetsColoringEachCategoryBuffer[iGroup]->getSize(),
	//			numThreads, &tetMeshGPU, stream);
	//}

}

void GAIA::PBDTetMeshNeoHookean::solveMaterialConstraintGPU_debugOnCPU(PBDTetMeshFEMGPU* pExternGPUTMeshOnCPU)
{
	// solveMaterialSequantialOnCPU((TetMeshFEMGPU_NeoHookean*)pExternGPUTMeshOnCPU, dt);

	for (size_t iParalGroup = 0; iParalGroup < tetsColoringEachCategoryBuffer.size(); iParalGroup++)
	{
		solveMaterialParallelOnCPU(tetsColoringEachCategoryBuffer[iParalGroup]->getCPUBuffer(), 
			tetsColoringEachCategoryBuffer[iParalGroup]->getSize(), (TetMeshFEMGPU_NeoHookean*)pExternGPUTMeshOnCPU);
	}
}

size_t GAIA::PBDTetMeshNeoHookean::numAllParallelizationGroups()
{
	return tetsColoringCategories().size();
}

size_t GAIA::PBDTetMeshNeoHookean::numTetParallelizationGroups()
{
	return PBDTetMeshFEMShared::numTetParallelizationGroups();
}

size_t GAIA::PBDTetMeshNeoHookean::getTetParallelizationGroupSize(int32_t iGroup)
{
	return PBDTetMeshFEMShared::getTetParallelizationGroupSize(iGroup);
}

int32_t* GAIA::PBDTetMeshNeoHookean::getTetParallelizationGroupBufferGPU(int32_t iGroup)
{
	return PBDTetMeshFEMShared::getTetParallelizationGroupBufferGPU(iGroup);
}



void GAIA::PBDTetMeshNeoHookean::initializeLambda()
{
	devLambdas = decltype(devLambdas)::Zero(devLambdas.size());
	volLambdas = decltype(volLambdas)::Zero(volLambdas.size());
}

GAIA::FloatingType GAIA::PBDTetMeshNeoHookean::evaluateElasticEnergy()
{
	FloatingType energy = 0.f;
	for (size_t iTet = 0; iTet < numTets(); iTet++)
	{
		int id0 = tetVIds()(0, iTet);
		int id1 = tetVIds()(1, iTet);
		int id2 = tetVIds()(2, iTet);
		int id3 = tetVIds()(3, iTet);

		Vec3 diff1 = vertex(id1) - vertex(id0);
		Vec3 diff2 = vertex(id2) - vertex(id0);
		Vec3 diff3 = vertex(id3) - vertex(id0);

		Eigen::Matrix<FloatingType, 3, 3> Dm;
		Dm << diff1[0], diff2[0], diff3[0],
			diff1[1], diff2[1], diff3[1],
			diff1[2], diff2[2], diff3[2]
			;

		int indToDsInv = iTet * 9;
		Eigen::Matrix<FloatingType, 3, 3> DmInv(DmInvs.data() + indToDsInv);
		Eigen::Matrix<FloatingType, 3, 3> F = Dm * DmInv;

		// Neo-Hookean Models
		FloatingType PHI_H = 0.f;
		if (pObjectParamsNeoHookean->volCompliance !=0)
		{
			PHI_H = (1 / pObjectParamsNeoHookean->volCompliance)* pow((F.determinant() - 1), 2);
		}
		FloatingType PHI_D = (1 / pObjectParamsNeoHookean->devCompliance) * ((F.transpose() * F).trace() - 3);
		energy += (PHI_H + PHI_D) * Dm.determinant() / 6;
	}
	
	return energy;
}



ObjectParamsPBDNeoHookean::SharedPtr GAIA::PBDTetMeshNeoHookean::getObjectParams()
{

	return std::static_pointer_cast<ObjectParamsPBDNeoHookean>(pObjectParams);
}

void GAIA::PBDTetMeshFEMShared_NeoHookean::copyAllToGPU()
{
	PBDTetMeshFEMShared::copyAllToGPU();
	devLambdasBuffer->toGPU();
	volLambdasBuffer->toGPU();
}


void GAIA::PBDTetMeshFEMShared_NeoHookean::initializeGPUTetMesh(PBDTetMeshFEM* pTetMesh, PBDTetMeshFEMGPU* pTetMeshGPU)
{
	PBDTetMeshFEMShared::initializeGPUTetMesh(pTetMesh, pTetMeshGPU);
	std::cout << "Initializing GPU tet mesh for NeoHookean material.\n";

	PBDTetMeshNeoHookean* pTetMeshPBDNeoHookean = (PBDTetMeshNeoHookean*)pTetMesh;
	devLambdasBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->numTets(), true, pTetMeshPBDNeoHookean->devLambdas.data());
	devLambdasBuffer->toGPU();
	volLambdasBuffer = std::make_shared<ManagedBuffer<FloatingTypeGPU>>(pTetMesh->numTets(), true, pTetMeshPBDNeoHookean->volLambdas.data());
	volLambdasBuffer->toGPU();

	TetMeshFEMGPU_NeoHookean* pTetMeshPBDNeoHookeanGPU = (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU;

	SET_GPU_MESH_BUFER_PTR(pTetMeshPBDNeoHookeanGPU, devLambdas);
	SET_GPU_MESH_BUFER_PTR(pTetMeshPBDNeoHookeanGPU, volLambdas);

	pTetMeshPBDNeoHookeanGPU->devCompliance = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->devCompliance;
	pTetMeshPBDNeoHookeanGPU->devDamping = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->devDamping;
	pTetMeshPBDNeoHookeanGPU->volCompliance = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->volCompliance;

}

void GAIA::PBDTetMeshFEMShared_NeoHookean::initializeGPUTetMeshForDebug(PBDTetMeshFEM* pTetMesh, PBDTetMeshFEMGPU* pTetMeshGPU)
{
	PBDTetMeshNeoHookean* pTetMeshPBDNeoHookean = (PBDTetMeshNeoHookean*)pTetMesh;
	TetMeshFEMGPU_NeoHookean* pTetMeshPBDNeoHookeanGPU = (TetMeshFEMGPU_NeoHookean*)pTetMeshGPU;

	PBDTetMeshFEMShared::initializeGPUTetMeshForDebug(pTetMesh, pTetMeshGPU);
	pTetMeshPBDNeoHookeanGPU->devCompliance = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->devCompliance;
	pTetMeshPBDNeoHookeanGPU->devDamping = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->devDamping;
	pTetMeshPBDNeoHookeanGPU->volCompliance = pTetMeshPBDNeoHookean->pObjectParamsNeoHookean->volCompliance;

	pTetMeshPBDNeoHookeanGPU->devLambdas = pTetMeshPBDNeoHookean->devLambdas.data();
	pTetMeshPBDNeoHookeanGPU->volLambdas = pTetMeshPBDNeoHookean->volLambdas.data();
}


//void GAIA::PBDTetMeshFEMShared_NeoHookean::createGPUTetMesh()
//{
//	std::cout << "Creating GPU tet mesh for NeoHookean material.\n";
//	pTetMeshGPU = std::make_shared<ManagedClassBuffer<TetMeshFEMGPU_NeoHookean>>();
//}


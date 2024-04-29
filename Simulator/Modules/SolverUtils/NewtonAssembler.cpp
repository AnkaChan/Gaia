#include "NewtonAssembler.h"

using namespace GAIA;

inline void GAIA::TetMeshNewtonAssembler::initialize(std::vector<TetMeshFEM::SharedPtr> meshes)
{
	numAllVertices = 0;
	numAllEdges = 0;
	numAllTets = 0;
	diagonalHessianBlockPtrs.resize(meshes.size());
	offDiagonalHessianBlockPtrs.resize(meshes.size());
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TetMeshFEM::SharedPtr pMesh = meshes[iMesh];

		numAllVertices += pMesh->numVertices();
		numAllEdges += pMesh->numEdges();
		numAllTets += pMesh->numTets();

		diagonalHessianBlockPtrs[iMesh].reserve(pMesh->numVertices() * 9);
		offDiagonalHessianBlockPtrs[iMesh].reserve(pMesh->numEdges() * 18);
	}
	newtonHessianTripletsElasticity.reserve(numAllVertices * 9 + numAllEdges * 18);
	int offset = 0;
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TetMeshFEM::SharedPtr pMesh = meshes[iMesh];

		//diagonal blocks
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPos + j, vertPos + i, 1.0);
				}
			}
		}

		//off-diagonal blocks
		for (int iE = 0; iE < pMesh->numEdges(); iE++)
		{
			int vertPosi = offset + pMesh->edges()(0, iE) * 3;
			int vertPosj = offset + pMesh->edges()(1, iE) * 3;
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPosi + i, vertPosj + j, 1.0);
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPosj + j, vertPosi + i, 1.0);
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}
	newtonHessianAll.resize(numAllVertices * 3, numAllVertices * 3);
	newtonHessianAll.setFromTriplets(newtonHessianTripletsElasticity.begin(), newtonHessianTripletsElasticity.end());
	newtonHessianAll.makeCompressed();

	if (solverType == 0)
	{
		solverDirect.analyzePattern(newtonHessianAll);
		
	}
	else {
		solverCG.compute(newtonHessianAll);
		solverCG.setMaxIterations(300);
		solverCG.setTolerance(1e-7f);
	}

	offset = 0;
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TetMeshFEM::SharedPtr pMesh = meshes[iMesh];
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					diagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianAll.coeffRef(vertPos + j, vertPos + i));
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
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianAll.coeffRef(vertPosi + i, vertPosj + j));
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianAll.coeffRef(vertPosj + j, vertPosi + i));
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}
	newtonForce.resize(numAllVertices * 3);
	elasticHessian.resize(numAllTets);
	elasticForce.resize(numAllTets);
}

void GAIA::TriMeshNewtonAssembler::initialize(std::vector<TriMeshFEM::SharedPtr> meshes_in, int solverType_in)
{
	meshes = meshes_in;
	solverType = solverType_in;

	numAllVertices = 0;
	numAllEdges = 0;
	numAllTris = 0;
	diagonalHessianBlockPtrs.resize(meshes.size());
	offDiagonalHessianBlockPtrs.resize(meshes.size());
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TriMeshFEM::SharedPtr pMesh = meshes[iMesh];

		numAllVertices += pMesh->numVertices();
		numAllEdges += pMesh->numEdges();
		numAllTris += pMesh->numFaces();

		diagonalHessianBlockPtrs[iMesh].reserve(pMesh->numVertices() * 9);
		offDiagonalHessianBlockPtrs[iMesh].reserve(pMesh->numEdges() * 18);
	}
	newtonHessianTripletsElasticity.reserve(numAllVertices * 9 + numAllEdges * 18);

	int offset = 0;
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TriMeshFEM::SharedPtr pMesh = meshes[iMesh];

		//diagonal blocks
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPos + j, vertPos + i, 1.0);
				}
			}
		}

		//off-diagonal blocks
		for (int iE = 0; iE < pMesh->numEdges(); iE++)
		{
			const EdgeInfo & edge = pMesh->getEdgeInfo(iE);
			int vertPosi = offset + edge.eV1 * 3;
			int vertPosj = offset + edge.eV2 * 3;
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPosi + i, vertPosj + j, 1.0);
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					newtonHessianTripletsElasticity.emplace_back(vertPosj + j, vertPosi + i, 1.0);
				}
			}
		}
		meshOffsets.push_back(offset);
		offset += pMesh->numVertices() * 3;
	}
	newtonHessianElasticity.resize(numAllVertices * 3, numAllVertices * 3);
	newtonHessianCollision.resize(numAllVertices * 3, numAllVertices * 3);
	newtonHessianAll.resize(numAllVertices * 3, numAllVertices * 3);

	newtonForce.resize(numAllVertices * 3);
	elasticHessian.resize(numAllTris);
	elasticForce.resize(numAllTris);
	elasticEnergy.resize(numAllTris);

	vfCollisionInfos.resize(meshes.size());
	eeCollisionInfos.resize(meshes.size());
	for (IdType iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		vfCollisionInfos[iMesh].resize(meshes[iMesh]->numVertices());
		eeCollisionInfos[iMesh].resize(meshes[iMesh]->numEdges());
	}

	newtonHessianElasticity.setFromTriplets(newtonHessianTripletsElasticity.begin(), newtonHessianTripletsElasticity.end());
	offset = 0;
	for (int iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		diagonalHessianBlockPtrs[iMesh].clear();
		offDiagonalHessianBlockPtrs[iMesh].clear();
		TriMeshFEM::SharedPtr pMesh = meshes[iMesh];
		for (int iV = 0; iV < pMesh->numVertices(); iV++)
		{
			int vertPos = offset + iV * 3;
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					diagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianElasticity.coeffRef(vertPos + j, vertPos + i));
				}
			}
		}
		for (int iE = 0; iE < pMesh->numEdges(); iE++)
		{
			const EdgeInfo& edge = pMesh->getEdgeInfo(iE);
			int vertPosi = offset + edge.eV1 * 3;
			int vertPosj = offset + edge.eV2 * 3;
			for (int j = 0; j < 3; ++j)
			{
				for (int i = 0; i < 3; ++i)
				{
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianElasticity.coeffRef(vertPosi + i, vertPosj + j));
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessianElasticity.coeffRef(vertPosj + j, vertPosi + i));
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}
}

void GAIA::TriMeshNewtonAssembler::analyzeCollision(const std::vector<std::vector<ClothVFContactQueryResult>>& vfCollisions, 
	const std::vector<std::vector<ClothEEContactQueryResult>>& eeCollisions)
{
	newtonHessianTripletsVFCollision.clear();
	newtonHessianTripletsEECollision.clear();

	size_t numSimulatedMeshes = meshes.size();
	for (IdType iMesh = 0; iMesh < vfCollisions.size(); iMesh++)
	{
		const std::vector<ClothVFContactQueryResult> vfCollisionsMesh = vfCollisions[iMesh];
		for (IdType iV = 0; iV < vfCollisionsMesh.size(); iV++)
		{
			const ClothVFContactQueryResult& vfCollision = vfCollisionsMesh[iV];
			for (IdType iVFContact = 0; iVFContact < vfCollision.contactPts.size(); iVFContact++)
			{
				const VFContactPointInfo& vfContact = vfCollision.contactPts[iVFContact];

				CIdType contactVId = vfContact.contactVertexId;
				CIdType contactVertexSideMeshId = vfContact.contactVertexSideMeshId;

				IdType t0 = -1, t1 = -1, t2 = -1;
				if (vfContact.contactFaceSideMeshId < numSimulatedMeshes)
				{
					const TriMeshFEM* pMeshFSide = meshes[vfContact.contactFaceSideMeshId].get();
					t0 = pMeshFSide->facePosVId(vfContact.contactFaceId, 0);
					t1 = pMeshFSide->facePosVId(vfContact.contactFaceId, 1);
					t2 = pMeshFSide->facePosVId(vfContact.contactFaceId, 2);
				}
				
				IdType vs[4] = { t0, t1, t2, contactVId };
				IdType meshIds[4] = {vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId, vfContact.contactFaceSideMeshId,contactVertexSideMeshId };

				for (IdType iRow = 0; iRow < 4; iRow++)
				{
					for (IdType iCol = 0; iCol < 4; iCol++)
					{
						if (vs[iCol] >=0 && vs[iRow] >= 0)
						{
							IdType vertPosRow = meshOffsets[meshIds[iRow]] + vs[iRow] * 3;
							IdType vertPosCol = meshOffsets[meshIds[iCol]] + vs[iCol] * 3;

							for (IdType i = 0; i < 3; ++i)
							{
								for (IdType j = 0; j < 3; ++j)
								{
									newtonHessianTripletsVFCollision.emplace_back(vertPosRow + i, vertPosCol + j, 1.0);
								}
							}
						}
					}
				}
			}
		}

		const std::vector<ClothEEContactQueryResult> eeCollisionsMesh = eeCollisions[iMesh];

		for (IdType iE = 0; iE < eeCollisionsMesh.size(); iE++)
		{
			const ClothEEContactQueryResult& eeCollision = eeCollisionsMesh[iE];
			for (IdType iEEContact = 0; iEEContact < eeCollision.contactPts.size(); iEEContact++)
			{
				const EEContactPointInfo& eeContact = eeCollision.contactPts[iEEContact];
				TriMeshFEM* pMesh = meshes[eeContact.contactMeshId1].get();


				IdType e2v1 = -1, e2v2 = -1;
				if (eeContact.contactMeshId2 < numSimulatedMeshes)
				{
					TriMeshFEM* pMeshOtherSide = meshes[eeContact.contactMeshId2].get();
					e2v1 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV1;
					e2v2 = pMeshOtherSide->getEdgeInfo(eeContact.contactEdgeId2).eV2;
				}

				IdType vs[4] = {
					pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV1, pMesh->getEdgeInfo(eeContact.contactEdgeId1).eV2,
					e2v1, e2v2
				};
				IdType meshIds[4] = { eeContact.contactMeshId1, eeContact.contactMeshId1, 
					eeContact.contactMeshId2, eeContact.contactMeshId2 };

				for (IdType iRow = 0; iRow < 2; iRow++)
				{
					// we only care about the half of the off-diagonal hessian namely hessian(this side, all side), 
					// the other half is handled by the collision results of other side
					for (IdType iCol = 0; iCol < 4; iCol++)
					{
						if (vs[iCol] >= 0 && vs[iRow] >= 0)
						{
							if (vs[iRow] >= 0 && vs[iCol]>=0)
							{
								IdType vertPosRow = meshOffsets[meshIds[iRow]] + vs[iRow] * 3;
								IdType vertPosCol = meshOffsets[meshIds[iCol]] + vs[iCol] * 3;

								for (IdType i = 0; i < 3; ++i)
								{
									for (IdType j = 0; j < 3; ++j)
									{
										newtonHessianTripletsEECollision.emplace_back(vertPosRow + i, vertPosCol + j, 1.0);
									}
								}
							}
						}
					}
				}
			}
		}
	}
	std::vector<NTriplet> newtonHessianTripletsAllCollisions = newtonHessianTripletsVFCollision;
	newtonHessianTripletsAllCollisions.insert(newtonHessianTripletsAllCollisions.end(), newtonHessianTripletsEECollision.begin(), newtonHessianTripletsEECollision.end());
	newtonHessianCollision.setFromTriplets(newtonHessianTripletsAllCollisions.begin(), newtonHessianTripletsAllCollisions.end());

}


void GAIA::TriMeshNewtonAssembler::updatePositions(const VecDynamic& dx)
{
	int offset = 0;
	for (size_t iMesh = 0; iMesh < meshes.size(); iMesh++)
	{
		TriMeshFEM* pMesh = meshes[iMesh].get();
		pMesh->positions() += Eigen::Map<const TVerticesMat>(dx.data() + offset, 3, pMesh->numVertices());
		offset += pMesh->numVertices() * 3;
	}
}


void GAIA::BaseNewtonAssembler::analyzePattern(bool makeCompressed)
{
	if (makeCompressed)
	{
		newtonHessianAll.makeCompressed();
	}

	if (solverType == 0)
	{
		solverDirect.analyzePattern(newtonHessianAll);
	}
	else if (solverType == 1)
	{
		solverCG.compute(newtonHessianAll);
		solverCG.setMaxIterations(cgMaxIterations);
		solverCG.setTolerance(cgTolerance);
	}
}

void GAIA::BaseNewtonAssembler::solve(bool patternChanged)
{
	newtonHessianAll = newtonHessianElasticity + newtonHessianCollision;
	if (patternChanged)
	{
#ifdef USE_MKL
		analyzePattern(true);
#else
		analyzePattern(false);
#endif // USE_MKL

	}

	if (solverType == 0)
	{
		solverDirect.factorize(newtonHessianAll);
		Ndx = solverDirect.solve(newtonForce);

		if (solverDirect.info() != Eigen::Success)
		{
			std::cerr << "Factorization failed. Error: " << solverDirect.info() << std::endl;
			std::exit(-1);
		}
	}
	else if (solverType == 1)
	{
		solverCG.compute(newtonHessianAll);
		Ndx = solverCG.solve(newtonForce);

		if (solverCG.info() != Eigen::Success)
		{
			std::cerr << "CG solve failed. Error: " << solverCG.info() << std::endl;
			//std::exit(-1);
		}
	} else
	{
		std::cout << "Error! Solver type not recognized!\n";
		exit(1);
	}

	if (Ndx.hasNaN()) {
		std::cerr << "Newton system has NaNs" << std::endl;
		std::exit(-1);
	}
}

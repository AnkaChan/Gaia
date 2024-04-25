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
	newtonHessian.resize(numAllVertices * 3, numAllVertices * 3);
	newtonHessian.setFromTriplets(newtonHessianTripletsElasticity.begin(), newtonHessianTripletsElasticity.end());
	newtonHessian.makeCompressed();

	if (solverType == 0)
	{
		solverDirect.analyzePattern(newtonHessian);
		
	}
	else {
		solverCG.compute(newtonHessian);
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
					diagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPos + j, vertPos + i));
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
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPosi + i, vertPosj + j));
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPosj + j, vertPosi + i));
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
		offset += pMesh->numVertices() * 3;
	}
	newtonHessian.resize(numAllVertices * 3, numAllVertices * 3);
	newtonForce.resize(numAllVertices * 3);
	elasticHessian.resize(numAllTris);
	elasticForce.resize(numAllTris);
	elasticEnergy.resize(numAllTris);
}

void GAIA::TriMeshNewtonAssembler::makeHessian(bool makeCompressed)
{
	

	std::vector<NTriplet> newtonHessianTripletsAll = newtonHessianTripletsElasticity;
	newtonHessianTripletsAll.insert(newtonHessianTripletsAll.end(), newtonHessianTripletsCollision.begin(), newtonHessianTripletsCollision.end());

	newtonHessian.setFromTriplets(newtonHessianTripletsElasticity.begin(), newtonHessianTripletsElasticity.end());
	if (makeCompressed)
	{
		newtonHessian.makeCompressed();
	}

	int offset = 0;
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
					diagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPos + j, vertPos + i));
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
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPosi + i, vertPosj + j));
				}
			}
			for (int i = 0; i < 3; ++i)
			{
				for (int j = 0; j < 3; ++j)
				{
					offDiagonalHessianBlockPtrs[iMesh].push_back(&newtonHessian.coeffRef(vertPosj + j, vertPosi + i));
				}
			}
		}
		offset += pMesh->numVertices() * 3;
	}

	if (solverType==0)
	{
		solverDirect.analyzePattern(newtonHessian);
	}
	else if (solverType == 1)
	{
		solverCG.compute(newtonHessian);
		solverCG.setMaxIterations(cgMaxIterations);
		solverCG.setTolerance(cgTolerance);
	}
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

void GAIA::BaseNewtonAssembler::solve()
{
	if (solverType == 0)
	{
		solverDirect.factorize(newtonHessian);
		Ndx = solverDirect.solve(newtonForce);

		if (solverDirect.info() != Eigen::Success)
		{
			std::cerr << "Factorization failed. Error: " << solverDirect.info() << std::endl;
			std::exit(-1);
		}
	}
	else if (solverType == 1)
	{
		solverCG.compute(newtonHessian);
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

#include "Test.h"
#include <common/math/vec2.h>
#include <common/math/vec3.h>
#include <common/math/vec4.h>

#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>

#include <VBDCloth/VBDTriMeshStVK.h>
#include <VBDCloth/VBDClothPhysics.h>

#include <Parallelization/CPUParallelization.h>

#include <MeshFrame/Memory/Array.h>

#include <set>


#include <CollisionDetector/TriMeshCollisionGeometry.h>
using namespace GAIA;

void testAtomicArray();
void Test_feasibleRegion();

#include <Parallelization/CPUParallelization.h>

void test()
{
	testAtomicArray();
	//Test_feasibleRegion();
}

void Test_feasibleRegion() {
	TriMeshFEM* pMesh = new TriMeshFEM();
	TriMeshParams::SharedPtr pObjParam = std::make_shared<TriMeshParams>();
	pObjParam->path = "F:/Code/02_Graphics/Energy-Based-Dynamics-Dev/CPP_Simulation/P14_APAPClothRegistration/TestData/Edge3.obj";
	pMesh->initialize(pObjParam, true);

	for (auto e : pMesh->pTopology->edgeInfos)
	{
		std::cout << "Edge: " << e.eV1 << ", " << e.eV2
			<< " | eV12Next: " << e.eV12Next
			<< " | eV21Next: " << e.eV21Next << std::endl;
	}

	for (size_t iNeiE = 0; iNeiE < 3; iNeiE++)
	{
		IdType edgeId = pMesh->pTopology->faces3NeighborEdges(iNeiE, 0);
		std::cout << "EdgeId: " << edgeId << "\n";
	}

	// edge feasible reagion
	embree::Vec3fa p1(-0.5f, 0.5f, -0.5f);

	bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);

	if (inEdgeFR )
	{
		std::cout << "p1 is in edge feasible region.\n";
	}

	for (int i = 0; i < 100; i++)
	{
		p1[2] += 0.01f ;
		if (i < 50)
		{
			bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
			std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";

		}
		else
		{
			bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
			std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";
		}
	}
	
	std::cout << "-----------------------------------------------\n";
	p1 = embree::Vec3fa(-0.5f, 0.5f, -0.5f);
	for (int i = 0; i < 100; i++)
	{
		p1[0] += 0.01f;
		bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
		std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";
	}

	std::cout << "-----------------------------------------------\n";
	p1 = embree::Vec3fa(-0.5f, 0.5f, -0.5f);
	for (int i = 0; i < 100; i++)
	{
		p1[0] += 0.01f;
		bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
		std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";
	}
	
	std::cout << "-----------------------------------------------\n";
	p1 = embree::Vec3fa(-0.5f, 0.5f, -0.5f);
	for (int i = 0; i < 100; i++)
	{
		p1[1] += 0.01f;
		bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
		std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";
	}

	std::cout << "-----------------------------------------------\n";
	p1 = embree::Vec3fa(0.5f, 0.5f, -0.5f);
	inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
	for (int i = 0; i < 100; i++)
	{
		p1[0] -= 0.01f;
		bool inEdgeFR = checkEdgeFeasibleRegion(p1, pMesh, 0, 2, 1e-7f);
		std::cout << "p1" << p1 << " is in edge feasible region: " << inEdgeFR << "\n";
	}

	// vertex feasible reagion
	//embree::Vec3fa p2(0, -0.5f, 0.f);

	//bool inVertexFR = checkVertexFeasibleRegion(p2, pMesh, 0, 1e-7f);
	//for (int i = 0; i < 100; i++)
	//{
	//	p2[1] += 0.01f;
	//	bool inVertFR = checkVertexFeasibleRegion(p2, pMesh, 0, 1e-7f);
	//	std::cout << "p2" << p2 << " is in vertex feasible region: " << inVertFR << "\n";
	//}

	//p2 = embree::Vec3fa(-0.1, -0.5f, 0.f);

	//std::cout << "-----------------------------------------------\n";
	//for (int i = 0; i < 100; i++)
	//{
	//	p2[0] += 0.01f;
	//	bool inVertFR = checkVertexFeasibleRegion(p2, pMesh, 0, 1e-7f);
	//	std::cout << "p2" << p2 << " is in vertex feasible region: " << inVertFR << "\n";
	//}
	//std::cout << "-----------------------------------------------\n";

	//p2 = embree::Vec3fa(0.1, 0.5f, 0.f);
	//p2 = embree::Vec3fa(0.0, 0.5f, 0.f);

	//std::cout << "-----------------------------------------------\n";
	//for (int i = 0; i < 100; i++)
	//{
	//	p2[1] += 0.01f;
	//	bool inVertFR = checkVertexFeasibleRegion(p2, pMesh, 2, 0);
	//	std::cout << "p2" << p2 << " is in vertex 2's feasible region: " << inVertFR << "\n";
	//}
	delete pMesh;

}

void testAtomicArray() {
	CPArrayStaticAtomic<int, 128> arr;
	int numThreads = 512;

	for (size_t i = 0; i < 10000; i++)
	{
		arr.clear();

		std::atomic<int> successCounter = 0;
		cpu_parallel_for(0, numThreads, [&](int i) {
			if (arr.push_back(i)) {
				successCounter ++;
			}
		});

		//std::cout << "arr.size(): " << arr.size() << " successCounter: " << successCounter << std::endl;
		//std::cout << "arr: " << arr << "\n";

		if (successCounter != arr.capacity() || arr.size() != arr.capacity())
		{
			std::cout << "arr.size(): " << arr.size() << " successCounter: " << successCounter << std::endl;
			std::cout << "arr: " << arr << "\n";
			getchar();
		}
	}
}

void testEdgesTopology(VBDTriMeshStVk::SharedPtr pStVKMesh)
{

	for (size_t iE = 0; iE < pStVKMesh->numEdges(); iE++)
	{
		std::cout << "Edge " << iE << ": " << pStVKMesh->pTopology->edgeInfos[iE].eV1 << " " << pStVKMesh->pTopology->edgeInfos[iE].eV2 << "\n";
		std::cout << "    eV12Next: " << pStVKMesh->pTopology->edgeInfos[iE].eV12Next << " | eV21Next: " << pStVKMesh->pTopology->edgeInfos[iE].eV21Next << "\n";
		std::cout << "    f1: " << pStVKMesh->pTopology->edgeInfos[iE].fId1 << " | f2: " << pStVKMesh->pTopology->edgeInfos[iE].fId2 << "\n";
		
		int fId1 = -1, fId2 = -1;
		for (size_t iF = 0; iF < pStVKMesh->numFaces(); iF++)
		{
			// check if the face matches the edge info
			std::set<int> fSet1 = {pStVKMesh->facePosVId(iF, 0), pStVKMesh->facePosVId(iF, 0) , pStVKMesh->facePosVId(iF, 2) };
		}
	}
}

void accumulateBendingAndForce(VBDTriMeshStVk::SharedPtr pStVKMesh, Eigen::MatrixXf& h, VecDynamic& force, bool skipHessian = false)
{
	const FloatingType bendingStiffness = pStVKMesh->objectParams().bendingStiffness;
	for (size_t iE = 0; iE < pStVKMesh->numEdges(); iE++)
	{
		const EdgeInfo& eInfo = pStVKMesh->pTopology->edgeInfos[iE];
		const Mat4& Q = pStVKMesh->edgeLaplacianQuadraticForms[iE];

		if (eInfo.fId2 == -1)
			// boundary edge
		{
			continue;
		}

		Eigen::Matrix<FloatingType, 4, 3> Xs;
		Xs.row(0) = pStVKMesh->vertex(eInfo.eV1);
		Xs.row(1) = pStVKMesh->vertex(eInfo.eV2);
		Xs.row(2) = pStVKMesh->vertex(eInfo.eV12Next);
		Xs.row(3) = pStVKMesh->vertex(eInfo.eV21Next);

		const Vec3 n1 = (Xs.row(1) - Xs.row(0)).cross(Xs.row(2) - Xs.row(0));
		const Vec3 n2 = (Xs.row(2) - Xs.row(3)).cross(Xs.row(1) - Xs.row(3));
		FloatingType n1_norm = n1.norm(), n2_norm = n2.norm();

		if (n1_norm < pStVKMesh->pPhysicsParams->degenerateTriangleThres
			|| n2_norm < pStVKMesh->pPhysicsParams->degenerateTriangleThres)
		{
			continue;
		}

		const Eigen::Matrix<FloatingType, 4, 3> dE_dXs = bendingStiffness * Q * Xs;

		VEC_BLOCK_DIM3x1_AT(force, eInfo.eV1) -= dE_dXs.row(0);
		VEC_BLOCK_DIM3x1_AT(force, eInfo.eV2) -= dE_dXs.row(1);
		VEC_BLOCK_DIM3x1_AT(force, eInfo.eV12Next) -= dE_dXs.row(2);
		VEC_BLOCK_DIM3x1_AT(force, eInfo.eV21Next) -= dE_dXs.row(3);

		if (!skipHessian)
		{
			int ids[4] = { eInfo.eV1, eInfo.eV2, eInfo.eV12Next, eInfo.eV21Next };
			for (size_t i = 0; i < 4; i++)
			{
				for (size_t j = 0; j < 4; j++) {
					int vId_i = ids[i];
					int vId_j = ids[j];
					h.block<3, 3>(3 * vId_i, 3 * vId_j) += bendingStiffness * Q(i, j) * Mat3::Identity();;
				}
			}

		}

	}
}

void accumulateStVKHessianAndForce(VBDTriMeshStVk::SharedPtr pStVKMesh, Eigen::MatrixXf & h, VecDynamic & force, bool skipHessian = false) 
{
	for (size_t iFace = 0; iFace < pStVKMesh->numFaces(); iFace++)
	{
		Mat3x2 F;
		pStVKMesh->calculateDeformationGradient(iFace, F);
		Mat2 G = pStVKMesh->greenStrain(F);

		CFloatingType lmbd = pStVKMesh->objectParams().lambda;
		CFloatingType miu = pStVKMesh->objectParams().miu;

		Mat2 S = 2.f * miu * G + lmbd * G.trace() * Mat2::Identity();

		auto DmInv = pStVKMesh->getDmInv(iFace);
		Mat3x2 F12 = -pStVKMesh->faceRestposeArea(iFace) * F * S * DmInv.transpose();

		int vId1 = pStVKMesh->facePosVId(iFace, 1);
		VEC_BLOCK_DIM3x1_AT(force, vId1) += F12.col(0);
		int vId2 = pStVKMesh->facePosVId(iFace, 2);
		VEC_BLOCK_DIM3x1_AT(force, vId2) += F12.col(1);
		int vId0 = pStVKMesh->facePosVId(iFace, 0);
		VEC_BLOCK_DIM3x1_AT(force, vId0) += -F12.col(0) - F12.col(1);

		if (skipHessian)
		{
			continue;
		}
		// compute Hessian diagonal blocks
		Eigen::Matrix<FloatingType, 6, 6> D2W_DFDF;

		const FloatingType Dm_inv1_1 = DmInv(0, 0);
		const FloatingType Dm_inv2_1 = DmInv(1, 0);
		const FloatingType Dm_inv1_2 = DmInv(0, 1);
		const FloatingType Dm_inv2_2 = DmInv(1, 1);

		const FloatingType F1_1 = F(0, 0);
		const FloatingType F2_1 = F(1, 0);
		const FloatingType F3_1 = F(2, 0);
		const FloatingType F1_2 = F(0, 1);
		const FloatingType F2_2 = F(1, 1);
		const FloatingType F3_2 = F(2, 1);

		const FloatingType F1_1_sqr = F1_1 * F1_1;
		const FloatingType F2_1_sqr = F2_1 * F2_1;
		const FloatingType F3_1_sqr = F3_1 * F3_1;
		const FloatingType F1_2_sqr = F1_2 * F1_2;
		const FloatingType F2_2_sqr = F2_2 * F2_2;
		const FloatingType F3_2_sqr = F3_2 * F3_2;

		const FloatingType e_uu = G(0, 0);
		const FloatingType e_vv = G(1, 1);
		const FloatingType e_uv = G(0, 1);
		const FloatingType e_uuvvSum = e_uu + e_vv;

		D2W_DFDF(0, 0) = F1_1 * (F1_1 * lmbd + 2 * F1_1 * miu) + 2 * miu * e_uu + lmbd * (
			e_uuvvSum)+F1_2_sqr * miu;

		D2W_DFDF(0, 1) = D2W_DFDF(1, 0) = F1_1 * (F2_1 * lmbd + 2 * F2_1 * miu) + F1_2 * F2_2 * miu;

		D2W_DFDF(0, 2) = D2W_DFDF(2, 0) = F1_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + F1_2 * F3_2 * miu;

		D2W_DFDF(0, 3) = D2W_DFDF(3, 0) = \
			2 * miu * e_uv + F1_1 * F1_2 * lmbd + F1_1 * F1_2 * miu;

		D2W_DFDF(0, 4) = D2W_DFDF(4, 0) = \
			F1_1 * F2_2 * lmbd + F1_2 * F2_1 * miu;

		D2W_DFDF(0, 5) = D2W_DFDF(5, 0) = \
			F1_1 * F3_2 * lmbd + F1_2 * F3_1 * miu;

		D2W_DFDF(1, 1) = \
			F2_1 * (F2_1 * lmbd + 2 * F2_1 * miu) + 2 * miu * e_uu + lmbd * (
				e_uuvvSum)+F2_2_sqr * miu;

		D2W_DFDF(1, 2) = D2W_DFDF(2, 1) = \
			F2_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + F2_2 * F3_2 * miu;

		D2W_DFDF(1, 3) = D2W_DFDF(3, 1) = \
			F1_2 * F2_1 * lmbd + F1_1 * F2_2 * miu;

		D2W_DFDF(1, 4) = D2W_DFDF(4, 1) = \
			2 * miu * e_uv + F2_1 * F2_2 * lmbd + F2_1 * F2_2 * miu;

		D2W_DFDF(1, 5) = D2W_DFDF(5, 1) = \
			F2_1 * F3_2 * lmbd + F2_2 * F3_1 * miu;

		D2W_DFDF(2, 2) = \
			F3_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + 2 * miu * e_uu + lmbd * (
				e_uuvvSum)+F3_2_sqr * miu;

		D2W_DFDF(2, 3) = D2W_DFDF(3, 2) = \
			F1_2 * F3_1 * lmbd + F1_1 * F3_2 * miu;

		D2W_DFDF(2, 4) = D2W_DFDF(4, 2) = \
			F2_2 * F3_1 * lmbd + F2_1 * F3_2 * miu;

		D2W_DFDF(2, 5) = D2W_DFDF(5, 2) = \
			2 * miu * e_uv + F3_1 * F3_2 * lmbd + F3_1 * F3_2 * miu;

		D2W_DFDF(3, 3) = \
			F1_2 * (F1_2 * lmbd + 2 * F1_2 * miu) + 2 * miu * e_vv + lmbd * (
				e_uuvvSum)+F1_1_sqr * miu;

		D2W_DFDF(3, 4) = D2W_DFDF(4, 3) = \
			F1_2 * (F2_2 * lmbd + 2 * F2_2 * miu) + F1_1 * F2_1 * miu;

		D2W_DFDF(3, 5) = D2W_DFDF(5, 3) = \
			F1_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + F1_1 * F3_1 * miu;

		D2W_DFDF(4, 4) = \
			F2_2 * (F2_2 * lmbd + 2 * F2_2 * miu) + 2 * miu * e_vv + lmbd * (
				e_uuvvSum)+F2_1_sqr * miu;

		D2W_DFDF(4, 5) = D2W_DFDF(5, 4) = \
			F2_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + F2_1 * F3_1 * miu;

		D2W_DFDF(5, 5) = \
			F3_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + 2 * miu * e_vv + lmbd * (
				e_uuvvSum)+F3_1_sqr * miu;

		D2W_DFDF = D2W_DFDF * pStVKMesh->faceRestposeArea(iFace);

		Eigen::Matrix < FloatingType, 6, 3> DF_DX1;
		DF_DX1 << -Dm_inv1_1 - Dm_inv2_1, 0, 0,
			0, -Dm_inv1_1 - Dm_inv2_1, 0,
			0, 0, -Dm_inv1_1 - Dm_inv2_1,
			-Dm_inv1_2 - Dm_inv2_2, 0, 0,
			0, -Dm_inv1_2 - Dm_inv2_2, 0,
			0, 0, -Dm_inv1_2 - Dm_inv2_2;

		Mat3 h_X1_X1 = DF_DX1.transpose() * D2W_DFDF * DF_DX1;



		Eigen::Matrix < FloatingType, 6, 3> DF_DX2;
		DF_DX2 << Dm_inv1_1, 0, 0,
			0, Dm_inv1_1, 0,
			0, 0, Dm_inv1_1,
			Dm_inv1_2, 0, 0,
			0, Dm_inv1_2, 0,
			0, 0, Dm_inv1_2;
		Mat3 h_X2_X2 = DF_DX2.transpose() * D2W_DFDF * DF_DX2;

		Eigen::Matrix < FloatingType, 6, 3> DF_DX3;
		DF_DX3 << Dm_inv2_1, 0, 0,
			0, Dm_inv2_1, 0,
			0, 0, Dm_inv2_1,
			Dm_inv2_2, 0, 0,
			0, Dm_inv2_2, 0,
			0, 0, Dm_inv2_2;
		//Mat3 h_X2_X3 = DF_DX3.transpose() * D2W_DFDF * DF_DX3;
		//std::cout << "D2W_DFDF:\n" << D2W_DFDF << "\n";

		std::array<Eigen::Matrix<FloatingType, 6, 3>, 3> DF_Di = { DF_DX1, DF_DX2, DF_DX3 };

		for (size_t i = 0; i < 3; i++)
		{
			int vId_i = pStVKMesh->facePosVId(iFace, i);
			for (size_t j = i; j < 3; j++) {
				int vId_j = pStVKMesh->facePosVId(iFace, j);

				Mat3 h_xi_xj = DF_Di[i].transpose() * D2W_DFDF * DF_Di[j];

				h.block<3, 3>(3 * vId_i, 3 * vId_j) += h_xi_xj;
				//std::cout << "h_x" <<i << "_h" << j <<":\n" << h_xi_xj << "\n";
				
				if (i!=j)
				{
					Mat3 h_xj_xi = DF_Di[j].transpose() * D2W_DFDF * DF_Di[i];
					h.block<3, 3>(3 * vId_j, 3 * vId_i) += h_xj_xi;
					//std::cout << "h_x" << j << "_h" << i << ":\n" << h_xj_xi << "\n";
				}
				// std::cout << "h:\n" << h;
				// std::cout << "\n";

			}
		}
	}

	// std::cout << "h:\n" << h << "\n";

}

void accumulateMaterialHessianAndForce(VBDTriMeshStVk::SharedPtr pStVKMesh, Eigen::MatrixXf& h, VecDynamic& force, bool skipHessian = false)
{
	accumulateStVKHessianAndForce(pStVKMesh, h, force, skipHessian);
	accumulateBendingAndForce(pStVKMesh, h, force, skipHessian);
}

void handleFixedPoints(const std::vector<int>& fixedPointList, Eigen::MatrixXf& h, VecDynamic& force) {
	for (size_t iFixedVert = 0; iFixedVert < fixedPointList.size(); iFixedVert++)
	{
		int vertId = fixedPointList[iFixedVert];
		VEC_BLOCK_DIM3x1_AT(force, vertId).setZero();

		h.row(vertId * 3).setZero();
		h.row(vertId * 3 + 1).setZero();
		h.row(vertId * 3 + 2).setZero();

		h.col(vertId * 3).setZero();
		h.col(vertId * 3 + 1).setZero();
		h.col(vertId * 3 + 2).setZero();
	}
}

void addInertiaHessianAndForce(VBDTriMeshStVk::SharedPtr pStVKMesh, FloatingType dt, Eigen::MatrixXf& h, VecDynamic & force) {
	FloatingType dtSqrReciprocal = 1.f / (dt * dt);
	for (size_t iV = 0; iV < pStVKMesh->numVertices(); iV++)
	{
		VEC_BLOCK_DIM3x1_AT(force, iV) -= pStVKMesh->vertexMass(iV) * (pStVKMesh->vertex(iV) - pStVKMesh->inertia.col(iV)) * (dtSqrReciprocal);
		h.block<3, 3>(3 * iV, 3 * iV) += Mat3::Identity() * pStVKMesh->vertexMass(iV) * (dtSqrReciprocal);
	} 
}

void forwardStep_intertia(VBDTriMeshStVk::SharedPtr pStVKMesh, FloatingType dt) {
	
	pStVKMesh->positionsPrev = pStVKMesh->positions();
	if (pStVKMesh->pPhysicsParams->associateGravityWithInertia)
	{
		pStVKMesh->velocities.colwise() += dt * pStVKMesh->pPhysicsParams->gravity;
	}
	pStVKMesh->handleVeclocityConstraint();
	//std::cout << "velocities:\n" << pStVKMesh->velocities << "\n";
	pStVKMesh->positions() += dt * pStVKMesh->velocities;
	//std::cout << "positions:\n" << pStVKMesh->positions() << "\n";
	pStVKMesh->inertia = pStVKMesh->positions();
}

FloatingType lineSearch(VBDTriMeshStVk::SharedPtr pStVKMesh, FloatingType E0, const VecDynamic& dx, FloatingType alpha, FloatingType c, FloatingType tau,
	int maxNumIters, FloatingType minStepSizeGD, FloatingType& best_meAll, FloatingType& best_meInertia, FloatingType& best_meElastic_stvk, FloatingType& best_meElastic_bending) {
	FloatingType m = 0.f;
	TVerticesMat orgPosistions = pStVKMesh->positions();

	m += dx.squaredNorm();

	FloatingType bestAlpha = 0.f;
	best_meAll = E0;

	ConvergenceStats statsTemp(1);

	for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
	{
		for (size_t iV = 0; iV < pStVKMesh->numVertices(); iV++)
		{
			pStVKMesh->positions().col(iV) = orgPosistions.col(iV) + alpha * VEC_BLOCK_DIM3x1_AT(dx, iV);
		}

		FloatingType meAll, meInertia, meElastic_stvk, meElastic_bending;
		meAll = pStVKMesh->evaluatedMeritEnergy(meInertia, meElastic_stvk, meElastic_bending);

		if (meAll < best_meAll)
		{
			bestAlpha = alpha;
			best_meAll = meAll;

			best_meInertia = meInertia;
			best_meElastic_stvk = meElastic_stvk;
			best_meElastic_bending = meElastic_bending;
		}

		// first Wolfe condition 
		if (meAll < E0 - alpha * c * m)
		{
			break;
		}
		else
		{
			alpha = alpha * tau;
		}
	}

	if (bestAlpha != alpha)
	{
		for (size_t iV = 0; iV < pStVKMesh->numVertices(); iV++)
		{
			pStVKMesh->positions().col(iV) = orgPosistions.col(iV) + bestAlpha * VEC_BLOCK_DIM3x1_AT(dx, iV);
		}
	}

	return bestAlpha;
}

void validateHessian(VBDTriMeshStVk::SharedPtr pStVKMesh, Eigen::MatrixXf& h) {
	VecDynamic orgForce;
	Eigen::MatrixXf h_null;
	orgForce.resize(3 * pStVKMesh->numVertices());
	orgForce.setZero();

	accumulateMaterialHessianAndForce(pStVKMesh, h_null, orgForce, true);

	TVerticesMat orgPos = pStVKMesh->positions();

	Eigen::MatrixXf h_numerical;
	h_numerical.resizeLike(h);
	h_numerical.setZero();

	VecDynamic force;
	force.resizeLike(orgForce);
	FloatingType delta = 1e-4f;

	for (size_t i = 0; i < 3 * pStVKMesh->numVertices(); i++)
	{
		force.setZero();
		int iV = i / 3;
		int iDim = i % 3;
		pStVKMesh->positions().col(iV)(iDim) += delta;

		accumulateMaterialHessianAndForce(pStVKMesh, h_null, force, true);

		h_numerical.col(i) = -(force - orgForce) / delta;

		//std::cout << "row " << i << " of h_numerical:\n" << h_numerical.col(i).transpose() << std::endl;
		//std::cout << "row " << i << " of h : \n" << h.col(i).transpose() << std::endl;


		FloatingType err = (h_numerical.col(i) - h.col(i)).norm();
		FloatingType relativeErr = err / h.col(i).norm();
		if (relativeErr > 0.01f) {
			std::cout << "Warning! Large relativeErr: " << relativeErr << "\n";

			// std::cout << "row " << i << " of h_numerical:\n" << h_numerical.col(i).transpose() << std::endl;
			// std::cout << "row " << i << " of h : \n" << h.col(i).transpose() << std::endl;
		}
		

		pStVKMesh->positions().col(iV) = orgPos.col(iV);

	}

}

;
FloatingType computeAvgForceNorm(VecDynamic& force) {
	FloatingType norm = 0.f;
	for (size_t iV = 0; iV < force.size() / 3; iV++)
	{
		norm += VEC_BLOCK_DIM3x1_AT(force, iV).norm();
	}

	return norm / (force.size() / 3);
}




void simulateClothMeshStVK()
{
	VBDObjectParamsTriMeshStVk::SharedPtr pObjParam = std::make_shared<VBDObjectParamsTriMeshStVk>();
	//pObjParam->path = "F:/Code/02_Graphics/Energy-Based-Dynamics-Dev/ParameterGen/Python_ParameterGen_APAPCloth/DataGen/Data/SyntheticData/C4.obj";
	pObjParam->path = "F:/Code/02_Graphics/Energy-Based-Dynamics-Dev/ParameterGen/Python_ParameterGen_APAPCloth/DataGen/Data/SyntheticData/C30.obj";
	//std::string outFolder = "E:/Data2/E_VBD_Results/00_ImplicitEuler";
	//std::string outFolder = "E:/Data2/E_VBD_Results/00_ImplicitEuler/withBending_0";
	std::string outFolder = "E:/Data2/E_VBD_Results/00_ImplicitEuler/C30_60FP_1Iter";
	std::string outFolderDebug = outFolder + "/Debug";

	MF::IO::createFolder(outFolderDebug);

	std::vector<int> fixedPointList = { 300, 329 };

	//pObjParam->path = "F:/Code/02_Graphics/Energy-Based-Dynamics-Dev/CPP_Simulation/P14_APAPClothRegistration/TestData/Edge1.obj";
	//std::vector<int> fixedPointList = { };

	//pObjParam->path = "F:/Code/02_Graphics/Energy-Based-Dynamics-Dev/CPP_Simulation/P15_VBDCloth/TestOneTriangle/OneTriangle.obj";
	//std::vector<int> fixedPointList = { };
	
	pObjParam->bendingStiffness = 100.f;
	pObjParam->miu = 1e6f;
	pObjParam->lambda = 1e6f;
	pObjParam->density = 0.02f;

	int numSteps = 300;

	pObjParam->fixedPoints = fixedPointList;

	VBDTriMeshStVk::SharedPtr pStVKMesh = std::make_shared<VBDTriMeshStVk>();

	VBDClothSimulationFramework physic;
	VBDClothPhysicsParameters::SharedPtr pPhysicParams = std::make_shared<VBDClothPhysicsParameters>();
	pPhysicParams->gravity << 0, -1000, 0;
	pPhysicParams->dt = 0.033333f;

	size_t numIter = 1;

	const FloatingType dt = 0.0166666f;

	physic.basePhysicsParams = pPhysicParams;

	pStVKMesh->initialize(pObjParam, &physic);
	std::cout << "vertexMass:\n" << pStVKMesh->vertexMass << "\n";

	// testEdgesTopology(pStVKMesh);
	Eigen::MatrixXf h(3*pStVKMesh->numVertices(), 3 * pStVKMesh->numVertices());
	h.setZero();

	VecDynamic force;
	force.resize(3 * pStVKMesh->numVertices());

	//std::vector<Vec3> deformationField(3);

	//deformationField[0] << 0.123, 0.645, 0.11;
	//deformationField[1] << 0.234, 0.534, 0.123;
	//deformationField[2] << 0.22, 0.156, 0.287;

	//for (size_t iV = 0; iV < 3; iV++)
	//{
	//	//pMesh->vertex(iV) += 0.1f * Vec3::Random();
	//	pStVKMesh->vertex(iV) += deformationField[iV];
	//}

	//accumulateStVKHessianAndForce(pStVKMesh, h, force);
	//std::cout << "h:\n" << h << "\n";
	//std::cout << "force:\n" << force << "\n";
	//validateHessian(pStVKMesh, h);

	//handleFixedPoints(fixedPointList, h, force);

	FloatingType meAll;
	FloatingType meInertia;
	FloatingType meElastic_stvk; 
	FloatingType meElastic_bending;

	for (size_t step = 0; step < numSteps; step++)
	{
		force.setZero();
		h.setZero();

		forwardStep_intertia(pStVKMesh, dt);
		
		pStVKMesh->saveAsPLY(outFolderDebug + "/A_step" + MF::STR::padToNew(std::to_string(step), 4, '0') 
			+ "_iter" + MF::STR::padToNew("0", 2, '0') + ".ply");

		for (size_t iIter = 0; iIter < numIter; iIter++)
		{
			force.setZero();
			h.setZero();

			accumulateMaterialHessianAndForce(pStVKMesh, h, force);

			// validateHessian(pStVKMesh, h);

			addInertiaHessianAndForce(pStVKMesh, dt, h, force);

			handleFixedPoints(fixedPointList, h, force);
			// std::cout << "h with intertia:\n" << h << "\n";

			VecDynamic dx = h.colPivHouseholderQr().solve(force);
			//std::cout << "dx:\n" << dx << "\n";

			FloatingType avgForceNorm = computeAvgForceNorm(force);

			FloatingType meAllInitial = pStVKMesh->evaluatedMeritEnergy(meInertia, meElastic_stvk, meElastic_bending);
			FloatingType stepSize = lineSearch(pStVKMesh, meAllInitial, dx, 1.0, 0.01f, 0.5f, 12, 0,
				meAll, meInertia, meElastic_stvk, meElastic_bending);

			force.setZero();
			accumulateMaterialHessianAndForce(pStVKMesh, Eigen::MatrixXf(), force, true);
			addInertiaHessianAndForce(pStVKMesh, dt, h, force);
			handleFixedPoints(fixedPointList, h, force);
			FloatingType avgForceNormNew = computeAvgForceNorm(force);


			std::cout << "step: " << step  << " iter: " << iIter << " | initial energy: " << meAllInitial
				<< " | energy after step: " << meAll 
				<< " | inertia: " << meInertia << " | elastic_stvk: " << meElastic_stvk << " | elastic_bending: " << meElastic_bending
				<< " | step size: " << stepSize 
				<< " | avg force norm before iter: " << avgForceNorm << " after iter: " << avgForceNormNew
				<< "\n";

			
			pStVKMesh->saveAsPLY(outFolderDebug + "/A_step" + MF::STR::padToNew(std::to_string(step), 4, '0')
				+ "_iter" + MF::STR::padToNew(std::to_string(iIter), 2, '0') + ".ply");
		}
		pStVKMesh->saveAsPLY(outFolder + "/A" + MF::STR::padToNew(std::to_string(step), 8, '0') + ".ply");

		pStVKMesh->velocities = (pStVKMesh->positions() - pStVKMesh->positionsPrev) / dt;
		
	}
	std::cout << "Done!\n";
}

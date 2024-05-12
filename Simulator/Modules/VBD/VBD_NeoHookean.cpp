#include "VBD_NeoHookean.h"
#include "../Parallelization/CPUParallelization.h"

using namespace GAIA;

// #define APPLY_LOCAL_LINE_SEARCH
// #define PSD_FILTERING

#undef PSD_FILTERING
void GAIA::VBDTetMeshNeoHookean::solverIteration()
{
	const size_t numVerts = numVertices();
	const size_t numColors = verticesColoringCategories().size();
	for (int iColor = 0; iColor < numColors; iColor++)
	{
		const std::vector<int32_t>& coloring = verticesColoringCategories()[iColor];
		const size_t sizeColoring = coloring.size();
		cpu_parallel_for(0, sizeColoring, [&](int iV) {
			int vId = coloring[iV];
			if (!fixedMask[vId])
			{
				VBDStep(vId);
			}
			});
	}

}

void GAIA::VBDTetMeshNeoHookean::evaluateInternalForce()
{
	cpu_parallel_for(0, numVertices(), [&](int iV) {
		Vec3 internalForce;
		internalForce.setZero();
		accumlateMaterialForce(iV, internalForce);
		vertexInternalForces.col(iV) = internalForce;
		});
}

void GAIA::VBDTetMeshNeoHookean::VBDStep(int iV)
{

	// assemble Hessian
	Mat3 h;
	Vec3 force;
	h.setZero();
	force.setZero();

	accumlateInertiaForceAndHessian(iV, force, h);
	accumlateMaterialForceAndHessian(iV, force, h);

	// add external force
	force += vertexExternalForces.col(iV);

	// line search
	if (force.squaredNorm() > CMP_EPSILON2)
	{
		// Solve linear system
		Vec3 descentDirection;
		FloatingType stepSize = pPhysicsParams->stepSize;
		FloatingType lineSearchShrinkFactor = pPhysicsParams->tau;
		bool solverSuccess = CuMatrix::solve3x3_psd_stable(h.data(), force.data(), descentDirection.data());

		if (!solverSuccess)
		{
			stepSize = pPhysicsParams->stepSizeGD;
			descentDirection = force;
			lineSearchShrinkFactor = pPhysicsParams->tau_GD;
			std::cout << "Solver failed at vertex " << iV << std::endl;
		}

		// line search
#ifdef APPLY_LOCAL_LINE_SEARCH
		FloatingType meInertia = 0;
		FloatingType meElastic_elastic = 0;
		FloatingType initialEnergy = evaluateVertexMeritEnergy(iV, meInertia, meElastic_elastic);
		FloatingType e = backTracingLineSearchVBD(iV, descentDirection, initialEnergy, stepSize, 0.f,
			lineSearchShrinkFactor, pPhysicsParams->backtracingLineSearchMaxIters);

		if (isnan(e))
		{
			assert(false);
		}
#else
		vertex(iV) += stepSize * descentDirection;
#endif // LOCAL_LINE_SEARCH

	}
}

FloatingType GAIA::VBDTetMeshNeoHookean::evaluateVertexMeritEnergy(int iV, FloatingType& meInertia, FloatingType& meElastic_elastic)
{

	// inertia
	meInertia = 0.5 * vertexMass(iV) * (vertex(iV) - inertia.col(iV)).squaredNorm() * pPhysicsParams->dtSqrReciprocal;

	// elastic
	meElastic_elastic = 0;
	const size_t numNeiTest = getNumVertexNeighborTets(iV);
	CFloatingType miu = ObjectParametersMaterial().miu;
	CFloatingType lmbd = ObjectParametersMaterial().lmbd;
	CFloatingType a = 1 + miu / lmbd;

	for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
	{
		IdType tetId = getVertexNeighborTet(iV, iNeiTet);

		CFloatingType A = tetRestVolume(tetId);
		auto DmInv = getDmInv(tetId);
		Mat3 Ds;
		computeDs(Ds, tetId);

		Mat3 F = Ds * DmInv;

		CFloatingType Phi_D = 0.5f * (F.cwiseProduct(F).sum() - 3.f);
		CFloatingType detF = F.determinant();
		CFloatingType Phi_H = 0.5f * SQR(detF - a);

		CFloatingType E = A * (miu * Phi_D + lmbd * Phi_H);

		meElastic_elastic += E;
	}

	return meInertia + meElastic_elastic;
}

inline void assembleVertexVForceAndHessian(const Vec9& dE_dF, const Mat9& d2E_dF, CFloatingType m1, CFloatingType m2, CFloatingType m3,
	Vec3& force, Mat3& h)
{
	CFloatingType A1 = dE_dF(0);
	CFloatingType A2 = dE_dF(1);
	CFloatingType A3 = dE_dF(2);
	CFloatingType A4 = dE_dF(3);
	CFloatingType A5 = dE_dF(4);
	CFloatingType A6 = dE_dF(5);
	CFloatingType A7 = dE_dF(6);
	CFloatingType A8 = dE_dF(7);
	CFloatingType A9 = dE_dF(8);

	force << A1 * m1 + A4 * m2 + A7 * m3,
		A2* m1 + A5 * m2 + A8 * m3,
		A3* m1 + A6 * m2 + A9 * m3;

	Eigen::Matrix<FloatingType, 3, 9> HL;

	HL.row(0) = d2E_dF.row(0) * m1 + d2E_dF.row(3) * m2 + d2E_dF.row(6) * m3;
	HL.row(1) = d2E_dF.row(1) * m1 + d2E_dF.row(4) * m2 + d2E_dF.row(7) * m3;
	HL.row(2) = d2E_dF.row(2) * m1 + d2E_dF.row(5) * m2 + d2E_dF.row(8) * m3;

	h.col(0) = HL.col(0) * m1 + HL.col(3) * m2 + HL.col(6) * m3;
	h.col(1) = HL.col(1) * m1 + HL.col(4) * m2 + HL.col(7) * m3;
	h.col(2) = HL.col(2) * m1 + HL.col(5) * m2 + HL.col(8) * m3;
}

template<typename T>
inline void assembleTetForceAndHessian(const Eigen::Vector<T, 9>& dE_dF, const Eigen::Matrix<T, 9, 9>& d2E_dF, const Eigen::Vector<T, 12>& m,
	Eigen::Vector<T, 12>& force, Eigen::Matrix<T, 12, 12>& h)
{
	Eigen::Map<const Eigen::Matrix3<T>> dedf(dE_dF.data());
	for (int i = 0; i < 4; ++i) {
		force.segment<3>(i * 3) = -(dedf * m.segment<3>(i * 3));
	}

	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j) {
			Eigen::Matrix<T, 3, 9> HL;

			const auto mi1 = m(3 * i);
			const auto mi2 = m(3 * i + 1);
			const auto mi3 = m(3 * i + 2);
			const auto mj1 = m(3 * j);
			const auto mj2 = m(3 * j + 1);
			const auto mj3 = m(3 * j + 2);
			HL.row(0) = d2E_dF.row(0) * mi1 + d2E_dF.row(3) * mi2 + d2E_dF.row(6) * mi3;
			HL.row(1) = d2E_dF.row(1) * mi1 + d2E_dF.row(4) * mi2 + d2E_dF.row(7) * mi3;
			HL.row(2) = d2E_dF.row(2) * mi1 + d2E_dF.row(5) * mi2 + d2E_dF.row(8) * mi3;

			h.block<3, 1>(3 * i, 3 * j) = HL.col(0) * mj1 + HL.col(3) * mj2 + HL.col(6) * mj3;
			h.block<3, 1>(3 * i, 3 * j + 1) = HL.col(1) * mj1 + HL.col(4) * mj2 + HL.col(7) * mj3;
			h.block<3, 1>(3 * i, 3 * j + 2) = HL.col(2) * mj1 + HL.col(5) * mj2 + HL.col(8) * mj3;
		}

	}

}


void GAIA::VBDTetMeshNeoHookean::accumlateMaterialForceAndHessian(int iV, Vec3& force, Mat3& hessian)
{
	const size_t numNeiTest = getNumVertexNeighborTets(iV);
	CFloatingType miu = ObjectParametersMaterial().miu;
	CFloatingType lmbd = ObjectParametersMaterial().lmbd;
	CFloatingType a = 1 + miu / lmbd;

	for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
	{
		IdType tetId = getVertexNeighborTet(iV, iNeiTet);

		CFloatingType A = tetRestVolume(tetId);

		Vec3 forceTet;
		Mat3 hessianTet;

		auto DmInv = getDmInv(tetId);

		Mat3 Ds;
		computeDs(Ds, tetId);
		//std::cout << "Ds:\n" << Ds << std::endl;

		Mat3 F = Ds * DmInv;
		//std::cout << "F:\n" << F << std::endl;

		CFloatingType detF = F.determinant();

		Eigen::Map<Vec9> dPhi_D_dF(F.data());

		// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

		CFloatingType F1_1 = F(0, 0);
		CFloatingType F2_1 = F(1, 0);
		CFloatingType F3_1 = F(2, 0);
		CFloatingType F1_2 = F(0, 1);
		CFloatingType F2_2 = F(1, 1);
		CFloatingType F3_2 = F(2, 1);
		CFloatingType F1_3 = F(0, 2);
		CFloatingType F2_3 = F(1, 2);
		CFloatingType F3_3 = F(2, 2);

		Vec9 ddetF_dF;
		ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
			F1_3* F3_2 - F1_2 * F3_3,
			F1_2* F2_3 - F1_3 * F2_2,
			F2_3* F3_1 - F2_1 * F3_3,
			F1_1* F3_3 - F1_3 * F3_1,
			F1_3* F2_1 - F1_1 * F2_3,
			F2_1* F3_2 - F2_2 * F3_1,
			F1_2* F3_1 - F1_1 * F3_2,
			F1_1* F2_2 - F1_2 * F2_1;

		// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

		Mat9 d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

		CFloatingType k = detF - a;
		d2E_dF_dF(0, 4) += k * F3_3;
		d2E_dF_dF(4, 0) += k * F3_3;
		d2E_dF_dF(0, 5) += k * -F2_3;
		d2E_dF_dF(5, 0) += k * -F2_3;
		d2E_dF_dF(0, 7) += k * -F3_2;
		d2E_dF_dF(7, 0) += k * -F3_2;
		d2E_dF_dF(0, 8) += k * F2_2;
		d2E_dF_dF(8, 0) += k * F2_2;

		d2E_dF_dF(1, 3) += k * -F3_3;
		d2E_dF_dF(3, 1) += k * -F3_3;
		d2E_dF_dF(1, 5) += k * F1_3;
		d2E_dF_dF(5, 1) += k * F1_3;
		d2E_dF_dF(1, 6) += k * F3_2;
		d2E_dF_dF(6, 1) += k * F3_2;
		d2E_dF_dF(1, 8) += k * -F1_2;
		d2E_dF_dF(8, 1) += k * -F1_2;

		d2E_dF_dF(2, 3) += k * F2_3;
		d2E_dF_dF(3, 2) += k * F2_3;
		d2E_dF_dF(2, 4) += k * -F1_3;
		d2E_dF_dF(4, 2) += k * -F1_3;
		d2E_dF_dF(2, 6) += k * -F2_2;
		d2E_dF_dF(6, 2) += k * -F2_2;
		d2E_dF_dF(2, 7) += k * F1_2;
		d2E_dF_dF(7, 2) += k * F1_2;

		d2E_dF_dF(3, 7) += k * F3_1;
		d2E_dF_dF(7, 3) += k * F3_1;
		d2E_dF_dF(3, 8) += k * -F2_1;
		d2E_dF_dF(8, 3) += k * -F2_1;

		d2E_dF_dF(4, 6) += k * -F3_1;
		d2E_dF_dF(6, 4) += k * -F3_1;
		d2E_dF_dF(4, 8) += k * F1_1;
		d2E_dF_dF(8, 4) += k * F1_1;

		d2E_dF_dF(5, 6) += k * F2_1;
		d2E_dF_dF(6, 5) += k * F2_1;
		d2E_dF_dF(5, 7) += k * -F1_1;
		d2E_dF_dF(7, 5) += k * -F1_1;

		d2E_dF_dF *= lmbd;

		d2E_dF_dF(0, 0) += miu;
		d2E_dF_dF(1, 1) += miu;
		d2E_dF_dF(2, 2) += miu;
		d2E_dF_dF(3, 3) += miu;
		d2E_dF_dF(4, 4) += miu;
		d2E_dF_dF(5, 5) += miu;
		d2E_dF_dF(6, 6) += miu;
		d2E_dF_dF(7, 7) += miu;
		d2E_dF_dF(8, 8) += miu;

		d2E_dF_dF *= A;

#ifdef PSD_FILTERING
		Eigen::JacobiSVD<Mat9> svd(d2E_dF_dF, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Vec9 eigenVals = svd.singularValues();

		Mat9 diagonalEV = Mat9::Zero();
		for (size_t iDim = 0; iDim < 9; iDim++)
		{
			if (eigenVals(iDim) > 0) diagonalEV(iDim, iDim) = eigenVals(iDim);
			else
			{
				//std::cout << "Negative Hessian ecountered: " << eigenVals.transpose() << "\n";
			}
			d2E_dF_dF = svd.matrixU() * diagonalEV * svd.matrixV().transpose();
		}
#endif // PSD_FILTERING



		Vec9 dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);



		CFloatingType DmInv1_1 = DmInv(0, 0);
		CFloatingType DmInv2_1 = DmInv(1, 0);
		CFloatingType DmInv3_1 = DmInv(2, 0);
		CFloatingType DmInv1_2 = DmInv(0, 1);
		CFloatingType DmInv2_2 = DmInv(1, 1);
		CFloatingType DmInv3_2 = DmInv(2, 1);
		CFloatingType DmInv1_3 = DmInv(0, 2);
		CFloatingType DmInv2_3 = DmInv(1, 2);
		CFloatingType DmInv3_3 = DmInv(2, 2);

		int vertedTetVId = getVertexNeighborTetVertexOrder(iV, iNeiTet);

		Eigen::Matrix<FloatingType, 9, 3> dF_dxi;
		Vec3 dE_dxi;
		Mat3 d2E_dxi_dxi;
		FloatingType m1, m2, m3;

		switch (vertedTetVId)
		{
		case 0:
			m1 = -DmInv1_1 - DmInv2_1 - DmInv3_1;
			m2 = -DmInv1_2 - DmInv2_2 - DmInv3_2;
			m3 = -DmInv1_3 - DmInv2_3 - DmInv3_3;
			break;

		case 1:
			m1 = DmInv1_1;
			m2 = DmInv1_2;
			m3 = DmInv1_3;
			break;
		case 2:
			m1 = DmInv2_1;
			m2 = DmInv2_2;
			m3 = DmInv2_3;
			break;

		case 3:
			m1 = DmInv3_1;
			m2 = DmInv3_2;
			m3 = DmInv3_3;
			break;
		default:
			break;
		}
		assembleVertexVForceAndHessian(dE_dF, d2E_dF_dF, m1, m2, m3, dE_dxi, d2E_dxi_dxi);

		force -= dE_dxi;
		hessian += d2E_dxi_dxi;

		//printf("dE_dF of tet %d:\n", tetId);
		//std::cout << dE_dF.transpose() << "\n";
		//printf("force of tet %d:\n-------------\n", tetId);
		//std::cout << force.transpose() << "\n";
	}
}

void GAIA::VBDTetMeshNeoHookean::accumlateMaterialForceAndHessian2(int iV, Vec3& force, Mat3& hessian)
{
	const size_t numNeiTest = getNumVertexNeighborTets(iV);
	const auto& material = ObjectParametersMaterial();
	CFloatingType miu = material.miu;
	CFloatingType lmbd = material.lmbd;
	CFloatingType a = 1 + miu / lmbd;
	Vec3 displacement = vertex(iV) - vertexPrevPos(iV);

	for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
	{
		IdType tetId = getVertexNeighborTet(iV, iNeiTet);

		CFloatingType A = tetRestVolume(tetId);

		Vec3 forceTet;
		Mat3 hessianTet;

		auto DmInv = getDmInv(tetId);

		Mat3 Ds;
		computeDs(Ds, tetId);
		//std::cout << "Ds:\n" << Ds << std::endl;

		Mat3 F = Ds * DmInv;
		//std::cout << "F:\n" << F << std::endl;

		CFloatingType detF = F.determinant();

		Eigen::Map<Vec9> dPhi_D_dF(F.data());

		// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

		CFloatingType F1_1 = F(0, 0);
		CFloatingType F2_1 = F(1, 0);
		CFloatingType F3_1 = F(2, 0);
		CFloatingType F1_2 = F(0, 1);
		CFloatingType F2_2 = F(1, 1);
		CFloatingType F3_2 = F(2, 1);
		CFloatingType F1_3 = F(0, 2);
		CFloatingType F2_3 = F(1, 2);
		CFloatingType F3_3 = F(2, 2);

		Vec9 ddetF_dF;
		ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
			F1_3* F3_2 - F1_2 * F3_3,
			F1_2* F2_3 - F1_3 * F2_2,
			F2_3* F3_1 - F2_1 * F3_3,
			F1_1* F3_3 - F1_3 * F3_1,
			F1_3* F2_1 - F1_1 * F2_3,
			F2_1* F3_2 - F2_2 * F3_1,
			F1_2* F3_1 - F1_1 * F3_2,
			F1_1* F2_2 - F1_2 * F2_1;

		// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

		Mat9 d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

		CFloatingType k = detF - a;
		d2E_dF_dF(0, 4) += k * F3_3;
		d2E_dF_dF(4, 0) += k * F3_3;
		d2E_dF_dF(0, 5) += k * -F2_3;
		d2E_dF_dF(5, 0) += k * -F2_3;
		d2E_dF_dF(0, 7) += k * -F3_2;
		d2E_dF_dF(7, 0) += k * -F3_2;
		d2E_dF_dF(0, 8) += k * F2_2;
		d2E_dF_dF(8, 0) += k * F2_2;

		d2E_dF_dF(1, 3) += k * -F3_3;
		d2E_dF_dF(3, 1) += k * -F3_3;
		d2E_dF_dF(1, 5) += k * F1_3;
		d2E_dF_dF(5, 1) += k * F1_3;
		d2E_dF_dF(1, 6) += k * F3_2;
		d2E_dF_dF(6, 1) += k * F3_2;
		d2E_dF_dF(1, 8) += k * -F1_2;
		d2E_dF_dF(8, 1) += k * -F1_2;

		d2E_dF_dF(2, 3) += k * F2_3;
		d2E_dF_dF(3, 2) += k * F2_3;
		d2E_dF_dF(2, 4) += k * -F1_3;
		d2E_dF_dF(4, 2) += k * -F1_3;
		d2E_dF_dF(2, 6) += k * -F2_2;
		d2E_dF_dF(6, 2) += k * -F2_2;
		d2E_dF_dF(2, 7) += k * F1_2;
		d2E_dF_dF(7, 2) += k * F1_2;

		d2E_dF_dF(3, 7) += k * F3_1;
		d2E_dF_dF(7, 3) += k * F3_1;
		d2E_dF_dF(3, 8) += k * -F2_1;
		d2E_dF_dF(8, 3) += k * -F2_1;

		d2E_dF_dF(4, 6) += k * -F3_1;
		d2E_dF_dF(6, 4) += k * -F3_1;
		d2E_dF_dF(4, 8) += k * F1_1;
		d2E_dF_dF(8, 4) += k * F1_1;

		d2E_dF_dF(5, 6) += k * F2_1;
		d2E_dF_dF(6, 5) += k * F2_1;
		d2E_dF_dF(5, 7) += k * -F1_1;
		d2E_dF_dF(7, 5) += k * -F1_1;

		d2E_dF_dF *= lmbd;

		//d2E_dF_dF(0, 0) += miu;
		//d2E_dF_dF(1, 1) += miu;
		//d2E_dF_dF(2, 2) += miu;
		//d2E_dF_dF(3, 3) += miu;
		//d2E_dF_dF(4, 4) += miu;
		//d2E_dF_dF(5, 5) += miu;
		//d2E_dF_dF(6, 6) += miu;
		//d2E_dF_dF(7, 7) += miu;
		//d2E_dF_dF(8, 8) += miu;

		d2E_dF_dF *= A;

#ifdef PSD_FILTERING
		Eigen::JacobiSVD<Mat9> svd(d2E_dF_dF, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Vec9 eigenVals = svd.singularValues();

		Mat9 diagonalEV = Mat9::Zero();
		for (size_t iDim = 0; iDim < 9; iDim++)
		{
			if (eigenVals(iDim) > 0) diagonalEV(iDim, iDim) = eigenVals(iDim);
			else
			{
				//std::cout << "Negative Hessian ecountered: " << eigenVals.transpose() << "\n";
			}
			d2E_dF_dF = svd.matrixU() * diagonalEV * svd.matrixV().transpose();
		}
#endif // PSD_FILTERING



		Vec9 dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);



		CFloatingType DmInv1_1 = DmInv(0, 0);
		CFloatingType DmInv2_1 = DmInv(1, 0);
		CFloatingType DmInv3_1 = DmInv(2, 0);
		CFloatingType DmInv1_2 = DmInv(0, 1);
		CFloatingType DmInv2_2 = DmInv(1, 1);
		CFloatingType DmInv3_2 = DmInv(2, 1);
		CFloatingType DmInv1_3 = DmInv(0, 2);
		CFloatingType DmInv2_3 = DmInv(1, 2);
		CFloatingType DmInv3_3 = DmInv(2, 2);

		int vertedTetVId = getVertexNeighborTetVertexOrder(iV, iNeiTet);

		Eigen::Matrix<FloatingType, 9, 3> dF_dxi;
		Vec3 dE_dxi;
		Mat3 d2E_dxi_dxi;
		FloatingType m1, m2, m3;

		switch (vertedTetVId)
		{
		case 0:
			m1 = -DmInv1_1 - DmInv2_1 - DmInv3_1;
			m2 = -DmInv1_2 - DmInv2_2 - DmInv3_2;
			m3 = -DmInv1_3 - DmInv2_3 - DmInv3_3;
			break;

		case 1:
			m1 = DmInv1_1;
			m2 = DmInv1_2;
			m3 = DmInv1_3;
			break;
		case 2:
			m1 = DmInv2_1;
			m2 = DmInv2_2;
			m3 = DmInv2_3;
			break;

		case 3:
			m1 = DmInv3_1;
			m2 = DmInv3_2;
			m3 = DmInv3_3;
			break;
		default:
			break;
		}
		assembleVertexVForceAndHessian(dE_dF, d2E_dF_dF, m1, m2, m3, dE_dxi, d2E_dxi_dxi);
		Mat3 dampingH = d2E_dxi_dxi * material.dampingHydrostatic;
		FloatingType tmp = (m1 * m1 + m2 * m2 + m3 * m3) * miu * A;
		d2E_dxi_dxi(0, 0) += tmp;
		d2E_dxi_dxi(1, 1) += tmp;
		d2E_dxi_dxi(2, 2) += tmp;
		tmp *= material.dampingDeviatoric;
		dampingH(0, 0) += tmp;
		dampingH(1, 1) += tmp;
		dampingH(2, 2) += tmp;
		dampingH /= pPhysicsParams->dt;
		Vec3 dampingForce = dampingH * displacement;
		force -= dE_dxi + dampingForce;
		hessian += d2E_dxi_dxi + dampingH;

		//printf("dE_dF of tet %d:\n", tetId);
		//std::cout << dE_dF.transpose() << "\n";
		//printf("force of tet %d:\n-------------\n", tetId);
		//std::cout << force.transpose() << "\n";
	}
}

// old version, use densit matmul
//void GAIA::VBDTetMeshNeoHookean::accumlateMaterialForceAndHessian(int iV, Vec3& force, Mat3& hessian)
//{
//	const size_t numNeiTest = getNumVertexNeighborTets(iV);
//	CFloatingType miu = ObjectParametersMaterial().miu;
//	CFloatingType lmbd = ObjectParametersMaterial().lmbd;
//	CFloatingType a = 1 + miu / lmbd;
//
//	for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
//	{
//		IdType tetId = getVertexNeighborTet(iV, iNeiTet);
//		
//		CFloatingType A = tetRestVolume(tetId);
//
//		Vec3 forceTet;
//		Mat3 hessianTet;
//
//		auto DmInv = getDmInv(tetId);
//
//		Mat3 Ds;
//		computeDs(Ds, tetId);
//
//		// std::cout << "Ds:\n" << Ds << std::endl;
//
//		Mat3 F = Ds * DmInv;
//		// std::cout << "F:\n" << F << std::endl;
//
//		// CFloatingType Phi_D = 0.5f * (F.cwiseProduct(F).sum() - 3.f);
//		CFloatingType detF = F.determinant();
//		// CFloatingType Phi_H = 0.5f * SQR(detF - a);
//
//		// CFloatingType E = A * (miu * Phi_D + lmbd * Phi_H);
//
//		// std::cout << "Phi_D: " << Phi_D << " | Phi_H: " << Phi_H << " | E: " << E << std::endl;
//
//		Eigen::Map<Vec9> dPhi_D_dF(F.data());
//
//		// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;
//
//		CFloatingType F1_1 = F(0, 0);
//		CFloatingType F2_1 = F(1, 0);
//		CFloatingType F3_1 = F(2, 0);
//		CFloatingType F1_2 = F(0, 1);
//		CFloatingType F2_2 = F(1, 1);
//		CFloatingType F3_2 = F(2, 1);
//		CFloatingType F1_3 = F(0, 2);
//		CFloatingType F2_3 = F(1, 2);
//		CFloatingType F3_3 = F(2, 2);
//
//		Vec9 ddetF_dF;
//		ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
//			F1_3* F3_2 - F1_2 * F3_3,
//			F1_2* F2_3 - F1_3 * F2_2,
//			F2_3* F3_1 - F2_1 * F3_3,
//			F1_1* F3_3 - F1_3 * F3_1,
//			F1_3* F2_1 - F1_1 * F2_3,
//			F2_1* F3_2 - F2_2 * F3_1,
//			F1_2* F3_1 - F1_1 * F3_2,
//			F1_1* F2_2 - F1_2 * F2_1;
//
//		// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;
//
//		Mat9 d2detF_dF_dF;
//		d2detF_dF_dF <<
//			0, 0, 0, 0, F3_3, -F2_3, 0, -F3_2, F2_2,
//			0, 0, 0, -F3_3, 0, F1_3, F3_2, 0, -F1_2,
//			0, 0, 0, F2_3, -F1_3, 0, -F2_2, F1_2, 0,
//			0, -F3_3, F2_3, 0, 0, 0, 0, F3_1, -F2_1,
//			F3_3, 0, -F1_3, 0, 0, 0, -F3_1, 0, F1_1,
//			-F2_3, F1_3, 0, 0, 0, 0, F2_1, -F1_1, 0,
//			0, F3_2, -F2_2, 0, -F3_1, F2_1, 0, 0, 0,
//			-F3_2, 0, F1_2, F3_1, 0, -F1_1, 0, 0, 0,
//			F2_2, -F1_2, 0, -F2_1, F1_1, 0, 0, 0, 0;
//
//		// std::cout << "d2detF_dF_dF:\n" << d2detF_dF_dF << std::endl;
//
//		Mat9 d2E_dF_dF = lmbd * (ddetF_dF * ddetF_dF.transpose() + (detF - a) * d2detF_dF_dF);
//
//		// Hessian of Phi_D
//		d2E_dF_dF(0, 0) += miu;
//		d2E_dF_dF(1, 1) += miu;
//		d2E_dF_dF(2, 2) += miu;
//		d2E_dF_dF(3, 3) += miu;
//		d2E_dF_dF(4, 4) += miu;
//		d2E_dF_dF(5, 5) += miu;
//		d2E_dF_dF(6, 6) += miu;
//		d2E_dF_dF(7, 7) += miu;
//		d2E_dF_dF(8, 8) += miu;
//
//		d2E_dF_dF *= A;
//
//		// std::cout << "d2E_dF_dF:\n" << d2E_dF_dF << std::endl;
//
//		Vec9 dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);
//
//		// std::cout << "dE_dF:\n" << dE_dF << std::endl;
//
//		CFloatingType DmInv1_1 = DmInv(0, 0);
//		CFloatingType DmInv2_1 = DmInv(1, 0);
//		CFloatingType DmInv3_1 = DmInv(2, 0);
//		CFloatingType DmInv1_2 = DmInv(0, 1);
//		CFloatingType DmInv2_2 = DmInv(1, 1);
//		CFloatingType DmInv3_2 = DmInv(2, 1);
//		CFloatingType DmInv1_3 = DmInv(0, 2);
//		CFloatingType DmInv2_3 = DmInv(1, 2);
//		CFloatingType DmInv3_3 = DmInv(2, 2);
//
//		int vertedTetVId = getVertexNeighborTetVertexOrder(iV, iNeiTet);
//
//		Eigen::Matrix<FloatingType, 9, 3> dF_dxi;
//		Vec3 dE_dxi;
//		Mat3 d2E_dxi_dxi;
//		switch (vertedTetVId)
//		{
//		case 0:
//			dF_dxi <<
//				-DmInv1_1 - DmInv2_1 - DmInv3_1, 0, 0,
//				0, -DmInv1_1 - DmInv2_1 - DmInv3_1, 0,
//				0, 0, -DmInv1_1 - DmInv2_1 - DmInv3_1,
//				-DmInv1_2 - DmInv2_2 - DmInv3_2, 0, 0,
//				0, -DmInv1_2 - DmInv2_2 - DmInv3_2, 0,
//				0, 0, -DmInv1_2 - DmInv2_2 - DmInv3_2,
//				-DmInv1_3 - DmInv2_3 - DmInv3_3, 0, 0,
//				0, -DmInv1_3 - DmInv2_3 - DmInv3_3, 0,
//				0, 0, -DmInv1_3 - DmInv2_3 - DmInv3_3;
//
//			dE_dxi = dE_dF.transpose() * dF_dxi;
//			d2E_dxi_dxi = dF_dxi.transpose() * d2E_dF_dF * dF_dxi;
//
//			// std::cout << "dE_dx1:\n" << dE_dxi << std::endl;
//			// std::cout << "d2E_dx1_dx1:\n" << d2E_dxi_dxi << std::endl;
//			break;
//
//		case 1:
//			dF_dxi <<
//				DmInv1_1, 0, 0,
//				0, DmInv1_1, 0,
//				0, 0, DmInv1_1,
//				DmInv1_2, 0, 0,
//				0, DmInv1_2, 0,
//				0, 0, DmInv1_2,
//				DmInv1_3, 0, 0,
//				0, DmInv1_3, 0,
//				0, 0, DmInv1_3;
//
//			dE_dxi = dE_dF.transpose() * dF_dxi;
//			d2E_dxi_dxi = dF_dxi.transpose() * d2E_dF_dF * dF_dxi;
//
//			// std::cout << "dE_dx2:\n" << dE_dxi << std::endl;
//			// std::cout << "d2E_dx2_dx2:\n" << d2E_dxi_dxi << std::endl;
//			break;
//		case 2:
//
//			dF_dxi <<
//				DmInv2_1, 0, 0,
//				0, DmInv2_1, 0,
//				0, 0, DmInv2_1,
//				DmInv2_2, 0, 0,
//				0, DmInv2_2, 0,
//				0, 0, DmInv2_2,
//				DmInv2_3, 0, 0,
//				0, DmInv2_3, 0,
//				0, 0, DmInv2_3;
//
//			dE_dxi = dE_dF.transpose() * dF_dxi;
//			d2E_dxi_dxi = dF_dxi.transpose() * d2E_dF_dF * dF_dxi;
//			
//			// std::cout << "dE_dx3:\n" << dF_dxi << std::endl;
//			// std::cout << "d2E_dx3_dx3:\n" << d2E_dxi_dxi << std::endl;
//			break;
//
//		case 3:
//			dF_dxi <<
//				DmInv3_1, 0, 0,
//				0, DmInv3_1, 0,
//				0, 0, DmInv3_1,
//				DmInv3_2, 0, 0,
//				0, DmInv3_2, 0,
//				0, 0, DmInv3_2,
//				DmInv3_3, 0, 0,
//				0, DmInv3_3, 0,
//				0, 0, DmInv3_3;
//
//			dE_dxi = dE_dF.transpose() * dF_dxi;
//			d2E_dxi_dxi = dF_dxi.transpose() * d2E_dF_dF * dF_dxi;
//
//			// std::cout << "dE_dx4:\n" << dE_dxi << std::endl;
//			// std::cout << "d2E_dx4_dx4:\n" << d2E_dxi_dxi << std::endl;
//			break;
//		default:
//			break;
//		}
//
//		force -= dE_dxi;
//		hessian += d2E_dxi_dxi;
//	}
//}

void GAIA::VBDTetMeshNeoHookean::accumlateMaterialForce(int iV, Vec3& force)
{
	const size_t numNeiTest = getNumVertexNeighborTets(iV);
	CFloatingType miu = ObjectParametersMaterial().miu;
	CFloatingType lmbd = ObjectParametersMaterial().lmbd;
	CFloatingType a = 1 + miu / lmbd;

	for (size_t iNeiTet = 0; iNeiTet < numNeiTest; iNeiTet++)
	{
		IdType tetId = getVertexNeighborTet(iV, iNeiTet);

		CFloatingType A = tetRestVolume(tetId);

		Vec3 forceTet;
		Mat3 hessianTet;

		auto DmInv = getDmInv(tetId);

		Mat3 Ds;
		computeDs(Ds, tetId);

		// std::cout << "Ds:\n" << Ds << std::endl;

		Mat3 F = Ds * DmInv;
		// std::cout << "F:\n" << F << std::endl;

		// CFloatingType Phi_D = 0.5f * (F.cwiseProduct(F).sum() - 3.f);
		CFloatingType detF = F.determinant();
		// CFloatingType Phi_H = 0.5f * SQR(detF - a);

		// CFloatingType E = A * (miu * Phi_D + lmbd * Phi_H);

		// std::cout << "Phi_D: " << Phi_D << " | Phi_H: " << Phi_H << " | E: " << E << std::endl;

		Eigen::Map<Vec9> dPhi_D_dF(F.data());

		// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

		CFloatingType F1_1 = F(0, 0);
		CFloatingType F2_1 = F(1, 0);
		CFloatingType F3_1 = F(2, 0);
		CFloatingType F1_2 = F(0, 1);
		CFloatingType F2_2 = F(1, 1);
		CFloatingType F3_2 = F(2, 1);
		CFloatingType F1_3 = F(0, 2);
		CFloatingType F2_3 = F(1, 2);
		CFloatingType F3_3 = F(2, 2);

		Vec9 ddetF_dF;
		ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
			F1_3* F3_2 - F1_2 * F3_3,
			F1_2* F2_3 - F1_3 * F2_2,
			F2_3* F3_1 - F2_1 * F3_3,
			F1_1* F3_3 - F1_3 * F3_1,
			F1_3* F2_1 - F1_1 * F2_3,
			F2_1* F3_2 - F2_2 * F3_1,
			F1_2* F3_1 - F1_1 * F3_2,
			F1_1* F2_2 - F1_2 * F2_1;

		Vec9 dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);

		// std::cout << "dE_dF:\n" << dE_dF << std::endl;

		CFloatingType DmInv1_1 = DmInv(0, 0);
		CFloatingType DmInv2_1 = DmInv(1, 0);
		CFloatingType DmInv3_1 = DmInv(2, 0);
		CFloatingType DmInv1_2 = DmInv(0, 1);
		CFloatingType DmInv2_2 = DmInv(1, 1);
		CFloatingType DmInv3_2 = DmInv(2, 1);
		CFloatingType DmInv1_3 = DmInv(0, 2);
		CFloatingType DmInv2_3 = DmInv(1, 2);
		CFloatingType DmInv3_3 = DmInv(2, 2);

		int vertedTetVId = getVertexNeighborTetVertexOrder(iV, iNeiTet);

		Eigen::Matrix<FloatingType, 9, 3> dF_dxi;
		Vec3 dE_dxi;
		switch (vertedTetVId)
		{
		case 0:
			dF_dxi <<
				-DmInv1_1 - DmInv2_1 - DmInv3_1, 0, 0,
				0, -DmInv1_1 - DmInv2_1 - DmInv3_1, 0,
				0, 0, -DmInv1_1 - DmInv2_1 - DmInv3_1,
				-DmInv1_2 - DmInv2_2 - DmInv3_2, 0, 0,
				0, -DmInv1_2 - DmInv2_2 - DmInv3_2, 0,
				0, 0, -DmInv1_2 - DmInv2_2 - DmInv3_2,
				-DmInv1_3 - DmInv2_3 - DmInv3_3, 0, 0,
				0, -DmInv1_3 - DmInv2_3 - DmInv3_3, 0,
				0, 0, -DmInv1_3 - DmInv2_3 - DmInv3_3;

			dE_dxi = dE_dF.transpose() * dF_dxi;

			// std::cout << "dE_dx1:\n" << dE_dxi << std::endl;
			// std::cout << "d2E_dx1_dx1:\n" << d2E_dxi_dxi << std::endl;
			break;

		case 1:
			dF_dxi <<
				DmInv1_1, 0, 0,
				0, DmInv1_1, 0,
				0, 0, DmInv1_1,
				DmInv1_2, 0, 0,
				0, DmInv1_2, 0,
				0, 0, DmInv1_2,
				DmInv1_3, 0, 0,
				0, DmInv1_3, 0,
				0, 0, DmInv1_3;

			dE_dxi = dE_dF.transpose() * dF_dxi;

			// std::cout << "dE_dx2:\n" << dE_dxi << std::endl;
			// std::cout << "d2E_dx2_dx2:\n" << d2E_dxi_dxi << std::endl;
			break;
		case 2:

			dF_dxi <<
				DmInv2_1, 0, 0,
				0, DmInv2_1, 0,
				0, 0, DmInv2_1,
				DmInv2_2, 0, 0,
				0, DmInv2_2, 0,
				0, 0, DmInv2_2,
				DmInv2_3, 0, 0,
				0, DmInv2_3, 0,
				0, 0, DmInv2_3;

			dE_dxi = dE_dF.transpose() * dF_dxi;

			// std::cout << "dE_dx3:\n" << dF_dxi << std::endl;
			// std::cout << "d2E_dx3_dx3:\n" << d2E_dxi_dxi << std::endl;
			break;

		case 3:
			dF_dxi <<
				DmInv3_1, 0, 0,
				0, DmInv3_1, 0,
				0, 0, DmInv3_1,
				DmInv3_2, 0, 0,
				0, DmInv3_2, 0,
				0, 0, DmInv3_2,
				DmInv3_3, 0, 0,
				0, DmInv3_3, 0,
				0, 0, DmInv3_3;

			dE_dxi = dE_dF.transpose() * dF_dxi;

			// std::cout << "dE_dx4:\n" << dE_dxi << std::endl;
			// std::cout << "d2E_dx4_dx4:\n" << d2E_dxi_dxi << std::endl;
			break;
		default:
			break;
		}

		force -= dE_dxi;
	}
}

void GAIA::VBDTetMeshNeoHookean::computeElasticForceHessian(int tetId, FloatingType& energy, Vec12& force, Mat12& hessian)
{
	const auto& material = ObjectParametersMaterial();
	CFloatingType miu = material.miu;
	CFloatingType lmbd = material.lmbd;
	CFloatingType a = 1 + miu / lmbd;

	CFloatingType A = tetRestVolume(tetId);

	auto DmInv = getDmInv(tetId);

	Mat3 Ds;
	computeDs(Ds, tetId);
	//std::cout << "Ds:\n" << Ds << std::endl;

	Mat3 F = Ds * DmInv;
	//std::cout << "F:\n" << F << std::endl;

	CFloatingType detF = F.determinant();

	Eigen::Map<Vec9> dPhi_D_dF(F.data());

	// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

	CFloatingType F1_1 = F(0, 0);
	CFloatingType F2_1 = F(1, 0);
	CFloatingType F3_1 = F(2, 0);
	CFloatingType F1_2 = F(0, 1);
	CFloatingType F2_2 = F(1, 1);
	CFloatingType F3_2 = F(2, 1);
	CFloatingType F1_3 = F(0, 2);
	CFloatingType F2_3 = F(1, 2);
	CFloatingType F3_3 = F(2, 2);

	Vec9 ddetF_dF;
	ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
		F1_3* F3_2 - F1_2 * F3_3,
		F1_2* F2_3 - F1_3 * F2_2,
		F2_3* F3_1 - F2_1 * F3_3,
		F1_1* F3_3 - F1_3 * F3_1,
		F1_3* F2_1 - F1_1 * F2_3,
		F2_1* F3_2 - F2_2 * F3_1,
		F1_2* F3_1 - F1_1 * F3_2,
		F1_1* F2_2 - F1_2 * F2_1;

	// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

	Mat9 d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

	CFloatingType k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	// d2E_dF_dF *= (lmbd * A);

	d2E_dF_dF *= lmbd;

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF *= A;

	energy = 0.5 * lmbd * k * k + 0.5 * miu * (F.squaredNorm() - 3);
	energy *= A;
	auto restEnergy = 0.5 * A * lmbd * SQR(1.f - a);
	energy -= restEnergy;

#ifdef PSD_FILTERING
	Eigen::JacobiSVD<Mat9> svd(d2E_dF_dF, Eigen::ComputeFullU | Eigen::ComputeFullV);

	Vec9 eigenVals = svd.singularValues();

	Mat9 diagonalEV = Mat9::Zero();
	for (size_t iDim = 0; iDim < 9; iDim++)
	{
		if (eigenVals(iDim) > 0) diagonalEV(iDim, iDim) = eigenVals(iDim);
		else
		{
			//std::cout << "Negative Hessian ecountered: " << eigenVals.transpose() << "\n";
		}
	}
	d2E_dF_dF = svd.matrixU() * diagonalEV * svd.matrixV().transpose();
#endif // PSD_FILTERING



	Vec9 dE_dF = A * (miu * dPhi_D_dF + (lmbd * detF - lmbd - miu) * ddetF_dF);



	CFloatingType DmInv1_1 = DmInv(0, 0);
	CFloatingType DmInv2_1 = DmInv(1, 0);
	CFloatingType DmInv3_1 = DmInv(2, 0);
	CFloatingType DmInv1_2 = DmInv(0, 1);
	CFloatingType DmInv2_2 = DmInv(1, 1);
	CFloatingType DmInv3_2 = DmInv(2, 1);
	CFloatingType DmInv1_3 = DmInv(0, 2);
	CFloatingType DmInv2_3 = DmInv(1, 2);
	CFloatingType DmInv3_3 = DmInv(2, 2);


	Vec12 m;

	m(0) = -DmInv1_1 - DmInv2_1 - DmInv3_1;
	m(1) = -DmInv1_2 - DmInv2_2 - DmInv3_2;
	m(2) = -DmInv1_3 - DmInv2_3 - DmInv3_3;
	m(3) = DmInv1_1;
	m(4) = DmInv1_2;
	m(5) = DmInv1_3;
	m(6) = DmInv2_1;
	m(7) = DmInv2_2;
	m(8) = DmInv2_3;
	m(9) = DmInv3_1;
	m(10) = DmInv3_2;
	m(11) = DmInv3_3;



	assembleTetForceAndHessian(dE_dF, d2E_dF_dF, m, force, hessian);
	Mat12 dampingH = hessian * material.dampingHydrostatic;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++)
	//	{
	//		FloatingType tmp = m.segment<3>(3 * i).dot(m.segment<3>(3 * j)) * miu * A;
	//		hessian(3 * i, 3 * j) += tmp;
	//		hessian(3 * i + 1, 3 * j + 1) += tmp;
	//		hessian(3 * i + 2, 3 * j + 2) += tmp;
	//		tmp *= material.dampingDeviatoric;
	//		dampingH(3 * i, 3 * j) += tmp;
	//		dampingH(3 * i + 1, 3 * j + 1) += tmp;
	//		dampingH(3 * i + 2, 3 * j + 2) += tmp;
	//	}
	//}
	dampingH /= pPhysicsParams->dt;
	Vec12 displacement;
	for (int i = 0; i < 4; ++i) {
		const auto iV = tetVIds()(i, tetId);
		displacement.segment<3>(3 * i) = vertex(iV) - vertexPrevPos(iV);
	}
	Vec12 dampingForce = dampingH * displacement;
	// energy += 0.5 * displacement.dot(dampingForce);
	force -= dampingForce;
	hessian += dampingH;

}

void GAIA::VBDTetMeshNeoHookean::computeElasticForceHessianDouble(int tetId, double& energy, Eigen::Vector<double, 12>& force, Eigen::Matrix<double, 12, 12>& hessian)
{
	const auto& material = ObjectParametersMaterial();
	const double miu = material.miu;
	const double lmbd = material.lmbd;
	const double a = 1 + miu / lmbd;

	const double A = tetRestVolume(tetId);

	Eigen::Matrix3d DmInv = getDmInv(tetId).cast<double>();

	Mat3 DsFloat;
	computeDs(DsFloat, tetId);
	Eigen::Matrix3d Ds = DsFloat.cast<double>();
	//std::cout << "Ds:\n" << Ds << std::endl;

	Eigen::Matrix3d F = Ds * DmInv;
	//std::cout << "F:\n" << F << std::endl;

	const double detF = F.determinant();

	Eigen::Map<const Eigen::Vector<double, 9>> dPhi_D_dF(F.data());

	// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

	const double F1_1 = F(0, 0);
	const double F2_1 = F(1, 0);
	const double F3_1 = F(2, 0);
	const double F1_2 = F(0, 1);
	const double F2_2 = F(1, 1);
	const double F3_2 = F(2, 1);
	const double F1_3 = F(0, 2);
	const double F2_3 = F(1, 2);
	const double F3_3 = F(2, 2);

	Eigen::Vector<double, 9> ddetF_dF;
	ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
		F1_3* F3_2 - F1_2 * F3_3,
		F1_2* F2_3 - F1_3 * F2_2,
		F2_3* F3_1 - F2_1 * F3_3,
		F1_1* F3_3 - F1_3 * F3_1,
		F1_3* F2_1 - F1_1 * F2_3,
		F2_1* F3_2 - F2_2 * F3_1,
		F1_2* F3_1 - F1_1 * F3_2,
		F1_1* F2_2 - F1_2 * F2_1;

	// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

	Eigen::Matrix<double, 9, 9> d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

	const double k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	// d2E_dF_dF *= (lmbd * A);

	d2E_dF_dF *= lmbd;

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF *= A;
	energy = 0.5 * lmbd * k * k + 0.5 * miu * (F.squaredNorm() - 3);
	energy *= A;
	auto restEnergy = 0.5 * A * lmbd * SQR(1.f - a);
	energy -= restEnergy;

#ifdef PSD_FILTERING
	Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(d2E_dF_dF, Eigen::ComputeFullU | Eigen::ComputeFullV);

	Eigen::Vector<double, 9> eigenVals = svd.singularValues();

	Eigen::Matrix<double, 9, 9> diagonalEV = Eigen::Matrix<double, 9, 9>::Zero();
	for (size_t iDim = 0; iDim < 9; iDim++)
	{
		if (eigenVals(iDim) > 0) diagonalEV(iDim, iDim) = eigenVals(iDim);
		else
		{
			//std::cout << "Negative Hessian ecountered: " << eigenVals.transpose() << "\n";
		}
	}
	d2E_dF_dF = svd.matrixU() * diagonalEV * svd.matrixV().transpose();
#endif // PSD_FILTERING



	Eigen::Vector<double, 9> dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);



	const auto DmInv1_1 = DmInv(0, 0);
	const auto DmInv2_1 = DmInv(1, 0);
	const auto DmInv3_1 = DmInv(2, 0);
	const auto DmInv1_2 = DmInv(0, 1);
	const auto DmInv2_2 = DmInv(1, 1);
	const auto DmInv3_2 = DmInv(2, 1);
	const auto DmInv1_3 = DmInv(0, 2);
	const auto DmInv2_3 = DmInv(1, 2);
	const auto DmInv3_3 = DmInv(2, 2);


	Eigen::Vector<double, 12> m;

	m(0) = -DmInv1_1 - DmInv2_1 - DmInv3_1;
	m(1) = -DmInv1_2 - DmInv2_2 - DmInv3_2;
	m(2) = -DmInv1_3 - DmInv2_3 - DmInv3_3;
	m(3) = DmInv1_1;
	m(4) = DmInv1_2;
	m(5) = DmInv1_3;
	m(6) = DmInv2_1;
	m(7) = DmInv2_2;
	m(8) = DmInv2_3;
	m(9) = DmInv3_1;
	m(10) = DmInv3_2;
	m(11) = DmInv3_3;



	assembleTetForceAndHessian(dE_dF, d2E_dF_dF, m, force, hessian);
	Eigen::Matrix<double, 12, 12> dampingH = hessian * material.dampingHydrostatic;
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++)
	//	{
	//		FloatingType tmp = m.segment<3>(3 * i).dot(m.segment<3>(3 * j)) * miu * A;
	//		hessian(3 * i, 3 * j) += tmp;
	//		hessian(3 * i + 1, 3 * j + 1) += tmp;
	//		hessian(3 * i + 2, 3 * j + 2) += tmp;
	//		tmp *= material.dampingDeviatoric;
	//		dampingH(3 * i, 3 * j) += tmp;
	//		dampingH(3 * i + 1, 3 * j + 1) += tmp;
	//		dampingH(3 * i + 2, 3 * j + 2) += tmp;
	//	}
	//}
	dampingH /= pPhysicsParams->dt;
	Eigen::Vector<double, 12> displacement;
	for (int i = 0; i < 4; ++i) {
		const auto iV = tetVIds()(i, tetId);
		displacement.segment<3>(3 * i) = (vertex(iV) - vertexPrevPos(iV)).cast<double>();
	}
	Eigen::Vector<double, 12> dampingForce = dampingH * displacement;
	// energy += 0.5 * displacement.dot(dampingForce);
	force -= dampingForce;
	hessian += dampingH;
}

void GAIA::VBDTetMeshNeoHookean::computeElasticGradientHessianF(const Mat3& F, FloatingType& energy, Vec9& gradient, Mat9& hessian)
{
	const auto& material = ObjectParametersMaterial();
	CFloatingType miu = material.miu;
	CFloatingType lmbd = material.lmbd;
	CFloatingType a = 1 + miu / lmbd;

	CFloatingType A = 1;

	//std::cout << "F:\n" << F << std::endl;

	CFloatingType detF = F.determinant();

	Eigen::Map<const Vec9> dPhi_D_dF(F.data());

	// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

	CFloatingType F1_1 = F(0, 0);
	CFloatingType F2_1 = F(1, 0);
	CFloatingType F3_1 = F(2, 0);
	CFloatingType F1_2 = F(0, 1);
	CFloatingType F2_2 = F(1, 1);
	CFloatingType F3_2 = F(2, 1);
	CFloatingType F1_3 = F(0, 2);
	CFloatingType F2_3 = F(1, 2);
	CFloatingType F3_3 = F(2, 2);

	Vec9 ddetF_dF;
	ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
		F1_3* F3_2 - F1_2 * F3_3,
		F1_2* F2_3 - F1_3 * F2_2,
		F2_3* F3_1 - F2_1 * F3_3,
		F1_1* F3_3 - F1_3 * F3_1,
		F1_3* F2_1 - F1_1 * F2_3,
		F2_1* F3_2 - F2_2 * F3_1,
		F1_2* F3_1 - F1_1 * F3_2,
		F1_1* F2_2 - F1_2 * F2_1;

	// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

	Mat9 d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

	CFloatingType k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	// d2E_dF_dF *= (lmbd * A);

	d2E_dF_dF *= lmbd;

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF *= A;
	energy = 0.5 * lmbd * k * k + 0.5 * miu * (F.squaredNorm() - 3);
	energy *= A;
	Vec9 dE_dF = A * (miu * dPhi_D_dF + (lmbd * detF - lmbd - miu) * ddetF_dF);
	gradient = dE_dF;
	hessian = d2E_dF_dF;
}

void GAIA::VBDTetMeshNeoHookean::computeElasticGradientHessianFDouble(const Eigen::Matrix3d& F, double& energy, Eigen::Vector<double, 9>& gradient, Eigen::Matrix<double, 9, 9>& hessian)
{
	const auto& material = ObjectParametersMaterial();
	const double miu = material.miu;
	const double lmbd = material.lmbd;
	const double a = 1 + miu / lmbd;

	const double A = 1;

	//std::cout << "F:\n" << F << std::endl;

	const double detF = F.determinant();

	Eigen::Map<const Eigen::Vector<double, 9>> dPhi_D_dF(F.data());

	// std::cout << "dPhi_D_dF:\n" << dPhi_D_dF << std::endl;

	const double F1_1 = F(0, 0);
	const double F2_1 = F(1, 0);
	const double F3_1 = F(2, 0);
	const double F1_2 = F(0, 1);
	const double F2_2 = F(1, 1);
	const double F3_2 = F(2, 1);
	const double F1_3 = F(0, 2);
	const double F2_3 = F(1, 2);
	const double F3_3 = F(2, 2);

	Eigen::Vector<double, 9> ddetF_dF;
	ddetF_dF << F2_2 * F3_3 - F2_3 * F3_2,
		F1_3* F3_2 - F1_2 * F3_3,
		F1_2* F2_3 - F1_3 * F2_2,
		F2_3* F3_1 - F2_1 * F3_3,
		F1_1* F3_3 - F1_3 * F3_1,
		F1_3* F2_1 - F1_1 * F2_3,
		F2_1* F3_2 - F2_2 * F3_1,
		F1_2* F3_1 - F1_1 * F3_2,
		F1_1* F2_2 - F1_2 * F2_1;

	// std::cout << "ddetF_dF:\n" << ddetF_dF << std::endl;

	Eigen::Matrix<double, 9, 9> d2E_dF_dF = ddetF_dF * ddetF_dF.transpose();

	const double k = detF - a;
	d2E_dF_dF(0, 4) += k * F3_3;
	d2E_dF_dF(4, 0) += k * F3_3;
	d2E_dF_dF(0, 5) += k * -F2_3;
	d2E_dF_dF(5, 0) += k * -F2_3;
	d2E_dF_dF(0, 7) += k * -F3_2;
	d2E_dF_dF(7, 0) += k * -F3_2;
	d2E_dF_dF(0, 8) += k * F2_2;
	d2E_dF_dF(8, 0) += k * F2_2;

	d2E_dF_dF(1, 3) += k * -F3_3;
	d2E_dF_dF(3, 1) += k * -F3_3;
	d2E_dF_dF(1, 5) += k * F1_3;
	d2E_dF_dF(5, 1) += k * F1_3;
	d2E_dF_dF(1, 6) += k * F3_2;
	d2E_dF_dF(6, 1) += k * F3_2;
	d2E_dF_dF(1, 8) += k * -F1_2;
	d2E_dF_dF(8, 1) += k * -F1_2;

	d2E_dF_dF(2, 3) += k * F2_3;
	d2E_dF_dF(3, 2) += k * F2_3;
	d2E_dF_dF(2, 4) += k * -F1_3;
	d2E_dF_dF(4, 2) += k * -F1_3;
	d2E_dF_dF(2, 6) += k * -F2_2;
	d2E_dF_dF(6, 2) += k * -F2_2;
	d2E_dF_dF(2, 7) += k * F1_2;
	d2E_dF_dF(7, 2) += k * F1_2;

	d2E_dF_dF(3, 7) += k * F3_1;
	d2E_dF_dF(7, 3) += k * F3_1;
	d2E_dF_dF(3, 8) += k * -F2_1;
	d2E_dF_dF(8, 3) += k * -F2_1;

	d2E_dF_dF(4, 6) += k * -F3_1;
	d2E_dF_dF(6, 4) += k * -F3_1;
	d2E_dF_dF(4, 8) += k * F1_1;
	d2E_dF_dF(8, 4) += k * F1_1;

	d2E_dF_dF(5, 6) += k * F2_1;
	d2E_dF_dF(6, 5) += k * F2_1;
	d2E_dF_dF(5, 7) += k * -F1_1;
	d2E_dF_dF(7, 5) += k * -F1_1;

	// d2E_dF_dF *= (lmbd * A);

	d2E_dF_dF *= lmbd;

	d2E_dF_dF(0, 0) += miu;
	d2E_dF_dF(1, 1) += miu;
	d2E_dF_dF(2, 2) += miu;
	d2E_dF_dF(3, 3) += miu;
	d2E_dF_dF(4, 4) += miu;
	d2E_dF_dF(5, 5) += miu;
	d2E_dF_dF(6, 6) += miu;
	d2E_dF_dF(7, 7) += miu;
	d2E_dF_dF(8, 8) += miu;

	d2E_dF_dF *= A;
	energy = 0.5 * lmbd * k * k + 0.5 * miu * (F.squaredNorm() - 3);
	energy *= A;
	Eigen::Vector<double, 9> dE_dF = A * (miu * dPhi_D_dF + lmbd * (detF - a) * ddetF_dF);
	gradient = dE_dF;
	hessian = d2E_dF_dF;
}

void GAIA::VBDTetMeshNeoHookean::validateElasticForceHessian(int tetId)
{
	double energy;
	Eigen::Vector<double, 12> force;
	Eigen::Matrix<double, 12, 12> hessian;
	Eigen::Vector<double, 12> forceNumeric;
	Eigen::Matrix<double, 12, 12> hessianNumeric;
	computeElasticForceHessianDouble(tetId, energy, force, hessian);
	const double eps = 1e-6;
	for (int i = 0; i < 4; ++i) {
		int vi = tetVIds()(i, tetId);
		for (int j = 0; j < 3; ++j) {
			FloatingType& x = mVertPos(j, vi);
			CFloatingType old_x = x;
			double real_eps = std::max<double>(1e-6f, eps * abs(x));
			x = double(x) + real_eps;
			double energyTmp;
			Eigen::Vector<double, 12> forceTmp;
			Eigen::Matrix<double, 12, 12> hessianTmp;
			computeElasticForceHessianDouble(tetId, energyTmp, forceTmp, hessianTmp);
			x = old_x;
			forceNumeric(i * 3 + j) = -(energyTmp - energy) / real_eps;
			hessianNumeric.col(i * 3 + j) = -(forceTmp - force) / real_eps;
		}
	}
	double forceDiff = (force - forceNumeric).norm() / 12;
	double hessianDiff = (hessian - hessianNumeric).norm() / 144;
	if (forceDiff > 1e-5) {
		std::cout << force.transpose() << std::endl;
		std::cout << forceNumeric.transpose() << std::endl;
	}
	if (hessianDiff > 1e-5) {
		std::cout << hessian << std::endl;
		std::cout << hessianNumeric << std::endl;
	}
}

void GAIA::VBDTetMeshNeoHookean::validateElasticGradientHessianF(int tetId)
{
	auto DmInv = getDmInv(tetId);
	Mat3 Ds;
	computeDs(Ds, tetId);
	Eigen::Matrix3d F = (Ds * DmInv).cast<double>();
	double energy;
	Eigen::Vector<double, 9> gradient;
	Eigen::Matrix<double, 9, 9> hessian;
	Eigen::Vector<double, 9> gradientNumeric;
	Eigen::Matrix<double, 9, 9> hessianNumeric;
	computeElasticGradientHessianFDouble(F, energy, gradient, hessian);
	const double eps = 1e-6;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			double& x = F(j, i);
			const double old_x = x;
			double real_eps = std::max<double>(1e-6, eps * abs(x));
			x = x + real_eps;
			double energyTmp;
			Eigen::Vector<double, 9> gradientTmp;
			Eigen::Matrix<double, 9, 9> hessianTmp;
			computeElasticGradientHessianFDouble(F, energyTmp, gradientTmp, hessianTmp);
			x = old_x;
			gradientNumeric(i * 3 + j) = (energyTmp - energy) / real_eps;
			hessianNumeric.col(i * 3 + j) = (gradientTmp - gradient) / real_eps;
		}
	}
	double gradientDiff = (gradient - gradientNumeric).norm() / 9;
	double hessianDiff = (hessian - hessianNumeric).norm() / 81;
	if (gradientDiff > 1e-5) {
		std::cout << gradient.transpose() << std::endl;
		std::cout << gradientNumeric.transpose() << std::endl;
	}
	if (hessianDiff > 1e-5) {
		std::cout << hessian << std::endl;
		std::cout << hessianNumeric << std::endl;
	}
}

bool GAIA::ObjectParametersVBDNeoHookean::fromJson(nlohmann::json& objectJsonParams)
{
	ObjectParamsVBD::fromJson(objectJsonParams);

	EXTRACT_FROM_JSON(objectJsonParams, miu);
	EXTRACT_FROM_JSON(objectJsonParams, lmbd);
	EXTRACT_FROM_JSON(objectJsonParams, initializationType);
	EXTRACT_FROM_JSON(objectJsonParams, dampingDeviatoric);
	EXTRACT_FROM_JSON(objectJsonParams, dampingHydrostatic);

	return true;
}

bool GAIA::ObjectParametersVBDNeoHookean::toJson(nlohmann::json& objectJsonParams)
{
	ObjectParamsVBD::toJson(objectJsonParams);
	PUT_TO_JSON(objectJsonParams, miu);
	PUT_TO_JSON(objectJsonParams, lmbd);
	PUT_TO_JSON(objectJsonParams, initializationType);
	PUT_TO_JSON(objectJsonParams, dampingDeviatoric);
	PUT_TO_JSON(objectJsonParams, dampingHydrostatic);

	return true;
}

void GAIA::VBDTetMeshNeoHookeanShared::initializeGPUTetMesh(VBDTetMeshNeoHookean* pTetMesh, VBDTetMeshNeoHookeanGPU* pTetMeshGPU, bool withAccelerator)
{
	VBDBaseTetMeshShared::initializeGPUTetMesh(pTetMesh, pTetMeshGPU, withAccelerator);

	pTetMeshGPU->miu = pTetMesh->ObjectParametersMaterial().miu;
	pTetMeshGPU->lmbd = pTetMesh->ObjectParametersMaterial().lmbd;
	pTetMeshGPU->dampingHydrostatic = pTetMesh->ObjectParametersMaterial().dampingHydrostatic;
	pTetMeshGPU->dampingDeviatoric = pTetMesh->ObjectParametersMaterial().dampingDeviatoric;

}

void GAIA::VBDTetMeshNeoHookeanShared::initializeGPUTetMeshForCPUDebug(TetMeshFEM* pTetMesh, TetMeshFEMGPU* pTetMeshGPU)
{
	VBDBaseTetMeshShared::initializeGPUTetMeshForCPUDebug(pTetMesh, pTetMeshGPU);

	VBDTetMeshNeoHookean* pTetMeshNeoHookean = (VBDTetMeshNeoHookean*)pTetMesh;
	VBDTetMeshNeoHookeanGPU* pTetMeshGPUNeoHookean = (VBDTetMeshNeoHookeanGPU*)pTetMeshGPU;

	pTetMeshGPUNeoHookean->miu = pTetMeshNeoHookean->ObjectParametersMaterial().miu;
	pTetMeshGPUNeoHookean->lmbd = pTetMeshNeoHookean->ObjectParametersMaterial().lmbd;

}

void VBDTetMeshNeoHookean::initializeGPUMesh()
{
	pTetMeshGPUBufferOnCPU = std::make_shared<VBDTetMeshNeoHookeanGPU>();

	pTetMeshShared = std::make_shared<VBDTetMeshNeoHookeanShared>();
	pTetMeshShared->initializeGPUTetMesh(this, pTetMeshGPUBufferOnCPU.get(), pPhysicsParams->useAccelerator);

	pTetMeshSharedBase = pTetMeshShared;

	pTetMeshGPUBuffer = std::make_shared<DeviceClassBuffer<VBDTetMeshNeoHookeanGPU>>();
	pTetMeshGPUBuffer->fromCPU(pTetMeshGPUBufferOnCPU.get());
}

VBDBaseTetMeshGPU* GAIA::VBDTetMeshNeoHookean::getGPUMesh()
{
	return pTetMeshGPUBuffer->getData();
}

VBDBaseTetMeshGPU* GAIA::VBDTetMeshNeoHookean::getGPUMesh_forCPUDebug()
{
	if (pTetMeshShared == nullptr)
	{
		pTetMeshShared = std::make_shared<VBDTetMeshNeoHookeanShared>();
		VBDTetMeshNeoHookeanGPU tetMeshGPUTempOnCPU;
		pTetMeshShared->initializeGPUTetMesh(this, &tetMeshGPUTempOnCPU, pPhysicsParams->useAccelerator);
		pTetMeshSharedBase = pTetMeshShared;
	}
	pTetMeshGPUBuffer_forCPUDebug = std::make_shared<VBDTetMeshNeoHookeanGPU>();
	pTetMeshShared->initializeGPUTetMeshForCPUDebug(this, pTetMeshGPUBuffer_forCPUDebug.get());

	return pTetMeshGPUBuffer_forCPUDebug.get();
}

void GAIA::VBDTetMeshNeoHookean::syncToGPU(bool sync, cudaStream_t stream)
{
	pTetMeshShared->syncToGPU(sync, stream);
}

void GAIA::VBDTetMeshNeoHookean::syncToCPU(bool sync, cudaStream_t stream)
{
	pTetMeshShared->syncToCPU(sync, stream);
}

void GAIA::VBDTetMeshNeoHookean::syncToGPUVertPosOnly(bool sync, cudaStream_t stream)
{
	pTetMeshShared->syncToGPUVertPosOnly(sync, stream);
}

void GAIA::VBDTetMeshNeoHookean::syncToCPUVertPosOnly(bool sync, cudaStream_t stream)
{
	pTetMeshShared->syncToCPUVertPosOnly(sync, stream);
}

void GAIA::VBDTetMeshNeoHookean::setGPUMeshActiveness(bool activeForCollision, bool activeForMaterialSolve)
{
	pTetMeshGPUBufferOnCPU->activeForCollision = activeForCollision;
	pTetMeshGPUBufferOnCPU->activeForMaterialSolve = activeForMaterialSolve;
	pTetMeshGPUBuffer->fromCPU(pTetMeshGPUBufferOnCPU.get());
}


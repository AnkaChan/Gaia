#include "VBDTriMeshStVK.h"
#include <Parallelization/CPUParallelization.h>

#include "VBDClothPhysicsParameters.h"

#define ARCCOS_THRES_HOLD (0.999f)

using namespace GAIA;

void GAIA::VBDTriMeshStVk::initialize(TriMeshParams::SharedPtr inObjectParams, BaseClothPhsicsFramework::Ptr in_pPhysics)
{
	VBDBaseTriMesh::initialize(inObjectParams, in_pPhysics);
	intializeEdgeQuadraticForms();
	targetPointsQueyResults.resize(numVertices());

}

//void GAIA::VBDTriMeshStVk::solve()
//{
//	//forwardStep_symplecticEuler();
//	//forwardStep_IE_GD();
//	forwardStep_IE_VBD();
//}



void GAIA::VBDTriMeshStVk::solverIteration()
{
	for (size_t iV = 0; iV < numVertices(); iV++)
	{
		VBDStep(iV);
	}
}

Vec3 GAIA::VBDTriMeshStVk::greenStrainVoigt(const Mat3x2& F)
{
	const Mat2 greenTensor = 0.5 * (F.transpose() * F - Mat2::Identity());
	return Vec3(greenTensor(0, 0), greenTensor(1, 1), 2 * greenTensor(0, 1));
}

Mat2 GAIA::VBDTriMeshStVk::greenStrain(const Mat3x2& F)
{
	return 0.5 * (F.transpose() * F - Mat2::Identity());
}

inline FloatingType getCotangent(const Vec3 v1, const Vec3 v2) {
	return v1.dot(v2) / v1.cross(v2).norm();
}

void GAIA::VBDTriMeshStVk::intializeEdgeQuadraticForms()
{
	edgeLaplacianQuadraticForms.resize(numEdges());

	cpu_parallel_for(0, numEdges(), [&](int iE) {
		const EdgeInfo& eInfo = pTopology->edgeInfos[iE];

		if (eInfo.fId2 == -1)
			// boundary edge
		{
			return;
		}

		Eigen::Matrix<FloatingType, 4, 3> Xs;
		Xs.row(0) = vertex(eInfo.eV1);
		Xs.row(1) = vertex(eInfo.eV2);
		Xs.row(2) = vertex(eInfo.eV12Next);
		Xs.row(3) = vertex(eInfo.eV21Next);

		const Vec3 e01 = Xs.row(1) - Xs.row(0); // e0
		const Vec3 e03 = Xs.row(3) - Xs.row(0); // e2
		const Vec3 e02 = Xs.row(2) - Xs.row(0); // e1
		const Vec3 e12 = Xs.row(2) - Xs.row(1); // e3
		const Vec3 e13 = Xs.row(3) - Xs.row(1); // e4

		CFloatingType c01 = getCotangent(e02, e01);
		CFloatingType c02 = getCotangent(e01, e03);
		CFloatingType c03 = getCotangent(-e01, e12);
		CFloatingType c04 = getCotangent(-e01, e13);

		Vec4 K0;
		K0 << c03 + c04, c01 + c02, -c01 - c03, -c02 - c04;

		CFloatingType A0 = computeNormal(eInfo.fId1, false).norm();
		CFloatingType A1 = computeNormal(eInfo.fId2, false).norm();

		edgeLaplacianQuadraticForms[iE] = 3.f * K0 * K0.transpose() / (A0 + A1);
		});
}

//void GAIA::VBDTriMeshStVk::evaluateInternalForce()
//{
//	evaluateInternalForcePerTriangle();
//	evaluateInternalForcePerEdge_quadratic();
//	evaluateInternalForcePerVertex();
//
//	handleInternalForceConstraint();
//
//}

//void GAIA::VBDTriMeshStVk::evaluateInternalForcePerTriangle()
//{
//	cpu_parallel_for(0, numFaces(), [&](int iFace) {
//		Mat3x2 F;
//		calculateDeformationGradient(iFace, F);
//		Mat2 G = greenStrain(F);
//
//		const CFloatingType lmbd = objectParams().lambda;
//		const CFloatingType miu = objectParams().miu;
//
//		const FloatingType faceArea = faceRestposeArea(iFace);
//
//		Mat2 S = 2.f * miu * G + lmbd * G.trace() * Mat2::Identity();
//
//		auto DmInv = getDmInv(iFace);
//		Mat3x2 F12 = -faceArea * F * S * DmInv.transpose();
//
//		TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 1, iFace) = F12.col(0);
//		TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 2, iFace) = F12.col(1);
//
//		TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 0, iFace) = -F12.col(0) - F12.col(1);
//
//		if (usePreconditioner)
//		{
//			// compute Hessian diagonal blocks
//			Eigen::Matrix<FloatingType, 6, 6> D2W_DFDF;
//
//			const FloatingType Dm_inv1_1 = DmInv(0, 0);
//			const FloatingType Dm_inv2_1 = DmInv(1, 0);
//			const FloatingType Dm_inv1_2 = DmInv(0, 1);
//			const FloatingType Dm_inv2_2 = DmInv(1, 1);
//
//			const FloatingType F1_1 = F(0, 0);
//			const FloatingType F2_1 = F(1, 0);
//			const FloatingType F3_1 = F(2, 0);
//			const FloatingType F1_2 = F(0, 1);
//			const FloatingType F2_2 = F(1, 1);
//			const FloatingType F3_2 = F(2, 1);
//
//			const FloatingType F1_1_sqr = F1_1 * F1_1;
//			const FloatingType F2_1_sqr = F2_1 * F2_1;
//			const FloatingType F3_1_sqr = F3_1 * F3_1;
//			const FloatingType F1_2_sqr = F1_2 * F1_2;
//			const FloatingType F2_2_sqr = F2_2 * F2_2;
//			const FloatingType F3_2_sqr = F3_2 * F3_2;
//
//			const FloatingType e_uu = G(0, 0);
//			const FloatingType e_vv = G(1, 1);
//			const FloatingType e_uv = G(0, 1);
//			const FloatingType e_uuvvSum = e_uu + e_vv;
//
//			D2W_DFDF(0, 0) = F1_1 * (F1_1 * lmbd + 2 * F1_1 * miu) + 2 * miu * e_uu + lmbd * (
//				e_uuvvSum)+F1_2_sqr * miu;
//
//			D2W_DFDF(0, 1) = D2W_DFDF(1, 0) = F1_1 * (F2_1 * lmbd + 2 * F2_1 * miu) + F1_2 * F2_2 * miu;
//
//			D2W_DFDF(0, 2) = D2W_DFDF(2, 0) = F1_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + F1_2 * F3_2 * miu;
//
//			D2W_DFDF(0, 3) = D2W_DFDF(3, 0) = \
//				2 * miu * e_uv + F1_1 * F1_2 * lmbd + F1_1 * F1_2 * miu;
//
//			D2W_DFDF(0, 4) = D2W_DFDF(4, 0) = \
//				F1_1 * F2_2 * lmbd + F1_2 * F2_1 * miu;
//
//			D2W_DFDF(0, 5) = D2W_DFDF(5, 0) = \
//				F1_1 * F3_2 * lmbd + F1_2 * F3_1 * miu;
//
//			D2W_DFDF(1, 1) = \
//				F2_1 * (F2_1 * lmbd + 2 * F2_1 * miu) + 2 * miu * e_uu + lmbd * (
//					e_uuvvSum)+F2_2_sqr * miu;
//
//			D2W_DFDF(1, 2) = D2W_DFDF(2, 1) = \
//				F2_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + F2_2 * F3_2 * miu;
//
//			D2W_DFDF(1, 3) = D2W_DFDF(3, 1) = \
//				F1_2 * F2_1 * lmbd + F1_1 * F2_2 * miu;
//
//			D2W_DFDF(1, 4) = D2W_DFDF(4, 1) = \
//				2 * miu * e_uv + F2_1 * F2_2 * lmbd + F2_1 * F2_2 * miu;
//
//			D2W_DFDF(1, 5) = D2W_DFDF(5, 1) = \
//				F2_1 * F3_2 * lmbd + F2_2 * F3_1 * miu;
//
//			D2W_DFDF(2, 2) = \
//				F3_1 * (F3_1 * lmbd + 2 * F3_1 * miu) + 2 * miu * e_uu + lmbd * (
//					e_uuvvSum)+F3_2_sqr * miu;
//
//			D2W_DFDF(2, 3) = D2W_DFDF(3, 2) = \
//				F1_2 * F3_1 * lmbd + F1_1 * F3_2 * miu;
//
//			D2W_DFDF(2, 4) = D2W_DFDF(4, 2) = \
//				F2_2 * F3_1 * lmbd + F2_1 * F3_2 * miu;
//
//			D2W_DFDF(2, 5) = D2W_DFDF(5, 2) = \
//				2 * miu * e_uv + F3_1 * F3_2 * lmbd + F3_1 * F3_2 * miu;
//
//			D2W_DFDF(3, 3) = \
//				F1_2 * (F1_2 * lmbd + 2 * F1_2 * miu) + 2 * miu * e_vv + lmbd * (
//					e_uuvvSum)+F1_1_sqr * miu;
//
//			D2W_DFDF(3, 4) = D2W_DFDF(4, 3) = \
//				F1_2 * (F2_2 * lmbd + 2 * F2_2 * miu) + F1_1 * F2_1 * miu;
//
//			D2W_DFDF(3, 5) = D2W_DFDF(5, 3) = \
//				F1_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + F1_1 * F3_1 * miu;
//
//			D2W_DFDF(4, 4) = \
//				F2_2 * (F2_2 * lmbd + 2 * F2_2 * miu) + 2 * miu * e_vv + lmbd * (
//					e_uuvvSum)+F2_1_sqr * miu;
//
//			D2W_DFDF(4, 5) = D2W_DFDF(5, 4) = \
//				F2_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + F2_1 * F3_1 * miu;
//
//			D2W_DFDF(5, 5) = \
//				F3_2 * (F3_2 * lmbd + 2 * F3_2 * miu) + 2 * miu * e_vv + lmbd * (
//					e_uuvvSum)+F3_1_sqr * miu;
//			D2W_DFDF = D2W_DFDF * faceArea;
//			Eigen::Matrix < FloatingType, 6, 3> DF_DX1;
//			DF_DX1 << -Dm_inv1_1 - Dm_inv2_1, 0, 0,
//				0, -Dm_inv1_1 - Dm_inv2_1, 0,
//				0, 0, -Dm_inv1_1 - Dm_inv2_1,
//				-Dm_inv1_2 - Dm_inv2_2, 0, 0,
//				0, -Dm_inv1_2 - Dm_inv2_2, 0,
//				0, 0, -Dm_inv1_2 - Dm_inv2_2;
//
//			triangleInternalForces_HessianDiagonalBlocks[iFace * 3] = DF_DX1.transpose() * D2W_DFDF * DF_DX1;
//
//			Eigen::Matrix < FloatingType, 6, 3> DF_DX2;
//			DF_DX2 << Dm_inv1_1, 0, 0,
//				0, Dm_inv1_1, 0,
//				0, 0, Dm_inv1_1,
//				Dm_inv1_2, 0, 0,
//				0, Dm_inv1_2, 0,
//				0, 0, Dm_inv1_2;
//			triangleInternalForces_HessianDiagonalBlocks[iFace * 3 + 1] = DF_DX2.transpose() * D2W_DFDF * DF_DX2;
//
//			Eigen::Matrix < FloatingType, 6, 3> DF_DX3;
//			DF_DX3 << Dm_inv2_1, 0, 0,
//				0, Dm_inv2_1, 0,
//				0, 0, Dm_inv2_1,
//				Dm_inv2_2, 0, 0,
//				0, Dm_inv2_2, 0,
//				0, 0, Dm_inv2_2;
//			triangleInternalForces_HessianDiagonalBlocks[iFace * 3 + 2] = DF_DX3.transpose() * D2W_DFDF * DF_DX3;
//		}
//
//		});
//}
//
//void GAIA::VBDTriMeshStVk::evaluateInternalForcePerVertex()
//{
//	cpu_parallel_for(0, numVertices(), [&](int iV) {
//		vertexInternalForces.col(iV).setZero();
//
//		for (size_t iNeiFace = 0; iNeiFace < numNeiFaces(iV); iNeiFace++)
//		{
//			int neiFaceId = getVertexIthNeiFace(iV, iNeiFace);
//			int neiFaceVertexOrder = getVertexIthNeiFaceOrder(iV, iNeiFace);
//
//			vertexInternalForces.col(iV) += TRIANGLE_INTERNAL_FORCE(triangleInternalForces, neiFaceVertexOrder, neiFaceId);
//			if (usePreconditioner)
//			{
//				vertexInternalForces_HessianDiagonalBlocks[iV] += triangleInternalForces_HessianDiagonalBlocks[neiFaceId * 3 + neiFaceVertexOrder];
//			}
//		}
//
//		});
//
//	// todo: parallelize this
//	for (size_t iE = 0; iE < pTopology->numEdges; iE++)
//	{
//		FloatingType bendingStiffness = objectParams().bendingStiffness;
//
//		const EdgeInfo& eInfo = pTopology->edgeInfos[iE];
//		const Mat4& Q = edgeLaplacianQuadraticForms[iE];
//
//		if (edgeDegenrateMask[iE])
//		{
//			continue;
//		}
//
//		vertexInternalForces.col(eInfo.eV1) += EDGE_INTERNAL_FORCE(edgeInternalForces, 0, iE);
//		vertexInternalForces.col(eInfo.eV2) += EDGE_INTERNAL_FORCE(edgeInternalForces, 1, iE);
//		vertexInternalForces.col(eInfo.eV12Next) += EDGE_INTERNAL_FORCE(edgeInternalForces, 2, iE);
//		vertexInternalForces.col(eInfo.eV21Next) += EDGE_INTERNAL_FORCE(edgeInternalForces, 3, iE);
//
//		if (usePreconditioner)
//		{
//			//vertexInternalForces_HessianDiagonalBlocks[eInfo.eV12Next] += edgeInternalForces_HessianDiagonalBlocks[iE * 4];
//			//vertexInternalForces_HessianDiagonalBlocks[eInfo.eV1] += edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 1];
//			//vertexInternalForces_HessianDiagonalBlocks[eInfo.eV2] += edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 2];
//			//vertexInternalForces_HessianDiagonalBlocks[eInfo.eV21Next] += edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 3];
//
//			vertexInternalForces_HessianDiagonalBlocks[eInfo.eV1] += bendingStiffness * Q(0, 0) * Mat3::Identity();;
//			vertexInternalForces_HessianDiagonalBlocks[eInfo.eV2] += bendingStiffness * Q(1, 1) * Mat3::Identity();;
//			vertexInternalForces_HessianDiagonalBlocks[eInfo.eV12Next] += bendingStiffness * Q(2, 2) * Mat3::Identity();
//			vertexInternalForces_HessianDiagonalBlocks[eInfo.eV21Next] += bendingStiffness * Q(3, 3) * Mat3::Identity();
//		}
//	}
//}

inline Mat3 mat3CrossProduct(const Mat3& M, const Vec3& a) {
	Mat3 result;
	result.col(0) = M.col(0).cross(a);
	result.col(1) = M.col(1).cross(a);
	result.col(2) = M.col(2).cross(a);

	return result;
}

//void GAIA::VBDTriMeshStVk::evaluateInternalForcePerEdge_bendingAngle()
//{
//	cpu_parallel_for(0, numEdges(), [&](int iE) {
//		const EdgeInfo& eInfo = pTopology->edgeInfos[iE];
//
//		if (eInfo.fId2 == -1)
//		{
//			edgeDegenrateMask[iE] = true;
//			edgeInternalForces.col(iE).setZero();
//			if (usePreconditioner)
//			{
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 1].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 2].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 3].setZero();
//			}
//			return;
//		}
//
//		const FloatingType bendingStiffness = objectParams().bendingStiffness;
//
//		const Vec3 x1 = vertex(eInfo.eV12Next);
//		const Vec3 x2 = vertex(eInfo.eV1);
//		const Vec3 x3 = vertex(eInfo.eV2);
//		const Vec3 x4 = vertex(eInfo.eV21Next);
//
//		const Vec3 n1 = (x2 - x1).cross(x3 - x1);
//		const Vec3 n2 = (x3 - x4).cross(x2 - x4);
//		FloatingType n1_norm = n1.norm(), n2_norm = n2.norm();
//
//		const Vec3 n1_n = n1 / n1_norm;
//		const Vec3 n2_n = n2 / n2_norm;
//
//		// degenerate triangle
//		// or flat bending
//		if (n1_norm < pPhysicsParams->degenerateTriangleThres
//			|| n2_norm < pPhysicsParams->degenerateTriangleThres)
//		{
//			edgeDegenrateMask[iE] = true;
//			edgeInternalForces.col(iE).setZero();
//			if (usePreconditioner)
//			{
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 1].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 2].setZero();
//				edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 3].setZero();
//			}
//			return;
//		}
//
//		const FloatingType cosTheta = n1_n.dot(n2_n);
//
//		const FloatingType sinHTheta_sqr = 0.5 * (1 - cosTheta);
//
//		const Vec3 dcosTheta_dx1 = 1 / n1_norm * (-(x3 - x1).cross(n2_n) + (x2 - x1).cross(n2_n));
//		const Vec3 dcosTheta_dx2 = 1 / n1_norm * (x3 - x1).cross(n2_n) - 1 / n2_norm * (x3 - x4).cross(n1_n);
//		const Vec3 dcosTheta_dx3 = -1 / n1_norm * (x2 - x1).cross(n2_n) + 1 / n2_norm * (x2 - x4).cross(n1_n);
//		const Vec3 dcosTheta_dx4 = 1 / n2_norm * (-(x2 - x4).cross(n1_n) + (x3 - x4).cross(n1_n));
//
//		const Vec3 dE_dx1 = -0.25 * bendingStiffness * dcosTheta_dx1;
//		const Vec3 dE_dx2 = -0.25 * bendingStiffness * dcosTheta_dx2;
//		const Vec3 dE_dx3 = -0.25 * bendingStiffness * dcosTheta_dx3;
//		const Vec3 dE_dx4 = -0.25 * bendingStiffness * dcosTheta_dx4;
//
//		edgeDegenrateMask[iE] = false;
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 0, iE) = -dE_dx1;
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 1, iE) = -dE_dx2;
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 2, iE) = -dE_dx3;
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 3, iE) = -dE_dx4;
//
//		if (usePreconditioner)
//		{
//			const Mat3 I3 = Mat3::Identity();
//
//			const Mat3 dn1_dx2 = mat3CrossProduct(I3, x3 - x1);
//			const Mat3 dn1_dx3 = -mat3CrossProduct(I3, x2 - x1);
//
//			const Mat3 dn2_dx2 = -mat3CrossProduct(I3, x3 - x4);
//			const Mat3 dn2_dx3 = mat3CrossProduct(I3, x2 - x4);
//
//			const Mat3 d2_cosTheta_dx_2_dx_2 = (1.f / (n1_norm * n2_norm)) * (mat3CrossProduct(dn1_dx2, x3 - x4) - mat3CrossProduct(dn2_dx2, x3 - x1));
//			const Mat3 d2_cosTheta_dx_3_dx_3 = (1.f / (n1_norm * n2_norm)) * (mat3CrossProduct(dn2_dx3, x2 - x1) - mat3CrossProduct(dn1_dx3, x2 - x4));
//
//			edgeInternalForces_HessianDiagonalBlocks[iE * 4].setZero(); // d2_E_dx_1_dx_1
//			edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 1] = -0.25 * bendingStiffness * d2_cosTheta_dx_2_dx_2;; // d2_E_dx_2_dx_2
//			edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 2] = -0.25 * bendingStiffness * d2_cosTheta_dx_3_dx_3; // d2_E_dx_3_dx_3
//			edgeInternalForces_HessianDiagonalBlocks[iE * 4 + 3].setZero(); // d2_E_dx_4_dx_4
//		}
//
//		// std::cout << "Force from edge " << iE << ": " << edgeInternalForces.col(iE).transpose() << std::endl;
//
//		});
////}

void GAIA::VBDTriMeshStVk::computeElasticityForceHessianForFace(int faceId, FloatingType& energy, Vec9& force, Mat9& hessian, bool psd)
{
	Mat3x2 F;
	calculateDeformationGradient(faceId, F);
	Mat2 G = greenStrain(F);

	const CFloatingType lmbd = objectParams().lambda;
	const CFloatingType miu = objectParams().miu;

	const FloatingType faceArea = faceRestposeArea(faceId);

	Mat2 S = 2.f * miu * G + lmbd * G.trace() * Mat2::Identity();

	auto DmInv = getDmInv(faceId);
	Mat3x2 F12 = -faceArea * F * S * DmInv.transpose();

	TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 1, faceId) = F12.col(0);
	TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 2, faceId) = F12.col(1);

	TRIANGLE_INTERNAL_FORCE(triangleInternalForces, 0, faceId) = -F12.col(0) - F12.col(1);

	// compute Hessian diagonal blocks
	Mat6 D2W_DFDF;

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
	D2W_DFDF = D2W_DFDF * faceArea;

	Eigen::Matrix < FloatingType, 6, 3> DF_DX1;
	DF_DX1 << -Dm_inv1_1 - Dm_inv2_1, 0, 0,
		0, -Dm_inv1_1 - Dm_inv2_1, 0,
		0, 0, -Dm_inv1_1 - Dm_inv2_1,
		-Dm_inv1_2 - Dm_inv2_2, 0, 0,
		0, -Dm_inv1_2 - Dm_inv2_2, 0,
		0, 0, -Dm_inv1_2 - Dm_inv2_2;


	Eigen::Matrix < FloatingType, 6, 3> DF_DX2;
	DF_DX2 << Dm_inv1_1, 0, 0,
		0, Dm_inv1_1, 0,
		0, 0, Dm_inv1_1,
		Dm_inv1_2, 0, 0,
		0, Dm_inv1_2, 0,
		0, 0, Dm_inv1_2;

	Eigen::Matrix < FloatingType, 6, 3> DF_DX3;
	DF_DX3 << Dm_inv2_1, 0, 0,
		0, Dm_inv2_1, 0,
		0, 0, Dm_inv2_1,
		Dm_inv2_2, 0, 0,
		0, Dm_inv2_2, 0,
		0, 0, Dm_inv2_2;

	Vec6 ms = {
		-Dm_inv1_1 - Dm_inv2_1,
		-Dm_inv1_2 - Dm_inv2_2,
		Dm_inv1_1,
		Dm_inv1_2,
		Dm_inv2_1,
		Dm_inv2_2
	};

	energy = faceArea * (0.5f * objectParams().lambda * pow(e_uu + e_vv, 2)
		+ objectParams().miu * (e_uu * e_uu + e_vv * e_vv + 2 * e_uv * e_uv));

	force.block<3, 1>(0, 0) = -F12.col(0) - F12.col(1);
	force.block<3, 1>(3, 0) = F12.col(0);
	force.block<3, 1>(6, 0) = F12.col(1);

	if (psd)
	{
		Eigen::JacobiSVD<Mat6> svd(D2W_DFDF, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Vec6 eigenVals = svd.singularValues();

		for (size_t iDim = 0; iDim < 6; iDim++)
		{
			if (eigenVals(iDim) < 0) eigenVals(iDim) = 0;
		}
		D2W_DFDF = svd.matrixU() * (eigenVals.asDiagonal()) * svd.matrixV().transpose();
	}
	assembleMembraneHessianForWholeFace(D2W_DFDF, ms, hessian);
}

void GAIA::VBDTriMeshStVk::computeElasticityEnergyForFace(int faceId, FloatingType& energy)
{
	Mat3x2 F;
	calculateDeformationGradient(faceId, F);
	Mat2 G = greenStrain(F);
	const FloatingType faceArea = faceRestposeArea(faceId);

	const FloatingType e_uu = G(0, 0);
	const FloatingType e_vv = G(1, 1);
	const FloatingType e_uv = G(0, 1);

	energy = faceArea * (0.5f * objectParams().lambda * pow(e_uu + e_vv, 2)
		+ objectParams().miu * (e_uu * e_uu + e_vv * e_vv + 2 * e_uv * e_uv));

}

//void GAIA::VBDTriMeshStVk::evaluateInternalForcePerEdge_quadratic()
//{
//	const FloatingType bendingStiffness = objectParams().bendingStiffness;
//
//	cpu_parallel_for(0, numEdges(), [&](int iE) {
//		const EdgeInfo& eInfo = pTopology->edgeInfos[iE];
//
//		if (eInfo.fId2 == -1)
//			// boundary edge
//		{
//			edgeDegenrateMask[iE] = true;
//			edgeInternalForces.col(iE).setZero();
//			return;
//		}
//
//		Eigen::Matrix<FloatingType, 4, 3> Xs;
//		Xs.row(0) = vertex(eInfo.eV1);
//		Xs.row(1) = vertex(eInfo.eV2);
//		Xs.row(2) = vertex(eInfo.eV12Next);
//		Xs.row(3) = vertex(eInfo.eV21Next);
//
//		const Vec3 n1 = (Xs.row(1) - Xs.row(0)).cross(Xs.row(2) - Xs.row(0));
//		const Vec3 n2 = (Xs.row(2) - Xs.row(3)).cross(Xs.row(1) - Xs.row(3));
//		FloatingType n1_norm = n1.norm(), n2_norm = n2.norm();
//
//		if (n1_norm < pPhysicsParams->degenerateTriangleThres
//			|| n2_norm < pPhysicsParams->degenerateTriangleThres)
//		{
//			edgeDegenrateMask[iE] = true;
//			edgeInternalForces.col(iE).setZero();
//			return;
//		}
//
//		const Eigen::Matrix<FloatingType, 4, 3> dE_dXs = bendingStiffness * edgeLaplacianQuadraticForms[iE] * Xs;
//
//		edgeDegenrateMask[iE] = false;
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 0, iE) = -dE_dXs.row(0);
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 1, iE) = -dE_dXs.row(1);
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 2, iE) = -dE_dXs.row(2);
//		EDGE_INTERNAL_FORCE(edgeInternalForces, 3, iE) = -dE_dXs.row(3);
//
//
//		});
//}

void GAIA::VBDTriMeshStVk::VBDStep(int iV)
{
	if (fixedMask(iV))
	{
		return;
	}
	// assemble Hessian
	Mat3 h;
	Vec3 force;
	h.setZero();
	force.setZero();

	accumlateInertiaForceAndHessian(iV, force, h);
	accumlateMaterialForceAndHessian(iV, force, h);

	// add external force
	force += vertexExternalForces.col(iV);

	// Solve linear system
	Vec3 descentDirection;
	FloatingType stepSize = pPhysicsParams->stepSize;
	FloatingType lineSearchShrinkFactor = 0.5f;
	bool solverSuccess = CuMatrix::solve3x3_psd_stable(h.data(), force.data(), descentDirection.data());

	if (!solverSuccess)
	{
		stepSize = pPhysicsParams->stepSizeGD;
		descentDirection = force;
		lineSearchShrinkFactor = 0.1f;
		std::cout << "Solver failed at vertex " << iV << std::endl;
		std::exit(-1);
	}

	// line search
	if (pPhysicsParams->useLineSearch && descentDirection.squaredNorm() > CMP_EPSILON2)
	{
		FloatingType meInertia = 0;
		FloatingType meElastic_stvk = 0;
		FloatingType meElastic_bending = 0;
		FloatingType initialEnergy = evaluateVertexMeritEnergy(iV, meInertia, meElastic_stvk, meElastic_bending);
		FloatingType e = backTracingLineSearchVBD(iV, descentDirection, initialEnergy, stepSize, 0.f, lineSearchShrinkFactor, pPhysicsParams->backtracingLineSearchMaxIters);

		if (isnan(e))
		{
			assert(false);
		}
	}
	else {
		vertex(iV) += stepSize * descentDirection;
	}
}

FloatingType GAIA::VBDTriMeshStVk::evaluateVertexMeritEnergy(int iV, FloatingType& meInertia, FloatingType& meElastic_stvk, FloatingType& meElastic_bending)
{
	meInertia = 0;
	meElastic_stvk = 0;
	meElastic_bending = 0;

	// inertia
	meInertia = 0.5 * vertexMass(iV) * (vertex(iV) - inertia.col(iV)).squaredNorm();

	// elastic
	for (size_t iNeiFace = 0; iNeiFace < numNeiFaces(iV); iNeiFace++)
	{
		int neiFaceId = getVertexIthNeiFace(iV, iNeiFace);
		Mat3x2 F;
		calculateDeformationGradient(neiFaceId, F);
		Mat2 G = greenStrain(F);
		FloatingType e_uu = G(0, 0), e_vv = G(1, 1), e_uv = G(1, 0);

		meElastic_stvk += faceRestposeArea(neiFaceId) * (0.5f * objectParams().lambda * pow(e_uu + e_vv, 2)
			+ objectParams().miu * (e_uu * e_uu + e_vv * e_vv + 2 * e_uv * e_uv));
	}

	meElastic_bending = 0.f;

	return meInertia + meElastic_stvk + meElastic_bending;
}

void GAIA::VBDTriMeshStVk::accumlateMaterialForceAndHessian(int iV, Vec3& grad, Mat3& hessian)
{
	accumlateStVKForceAndHessian(iV, grad, hessian);
	accumlateBendingForceAndHessian(iV, grad, hessian);
}

void GAIA::VBDTriMeshStVk::accumlateStVKForceAndHessian(int iV, Vec3& force, Mat3& hessian)
{
	size_t numberNeiFaces = numNeiFaces(iV);
	Mat3 old_hessian = hessian;

	for (size_t iNeiF = 0; iNeiF < numberNeiFaces; iNeiF++)
	{
		int faceId = getVertexIthNeiFace(iV, iNeiF);
		IdType faceVertexOrder = getVertexIthNeiFaceOrder(iV, iNeiF);

		Mat3x2 F;
		calculateDeformationGradient(faceId, F);
		Mat2 G = greenStrain(F);

		const CFloatingType lmbd = objectParams().lambda;
		const CFloatingType miu = objectParams().miu;

		const FloatingType faceArea = faceRestposeArea(faceId);

		Mat2 S = 2.f * miu * G + lmbd * G.trace() * Mat2::Identity();

		auto DmInv = getDmInv(faceId);
		Mat3x2 F12 = -faceArea * F * S * DmInv.transpose();

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
		D2W_DFDF = D2W_DFDF * faceArea;

		//if (iV == facePosVId(faceId, 0)) //Face vertex 1
		//{
		//	force += -F12.col(0) - F12.col(1);
		//	Eigen::Matrix < FloatingType, 6, 3> DF_DX1;
		//	DF_DX1 << -Dm_inv1_1 - Dm_inv2_1, 0, 0,
		//		0, -Dm_inv1_1 - Dm_inv2_1, 0,
		//		0, 0, -Dm_inv1_1 - Dm_inv2_1,
		//		-Dm_inv1_2 - Dm_inv2_2, 0, 0,
		//		0, -Dm_inv1_2 - Dm_inv2_2, 0,
		//		0, 0, -Dm_inv1_2 - Dm_inv2_2;

		//	hessian += DF_DX1.transpose() * D2W_DFDF * DF_DX1;
		//}
		//else if (iV == facePosVId(faceId, 1)) //Face vertex 2
		//{
		//	force += F12.col(0);

		//	Eigen::Matrix < FloatingType, 6, 3> DF_DX2;
		//	DF_DX2 << Dm_inv1_1, 0, 0,
		//		0, Dm_inv1_1, 0,
		//		0, 0, Dm_inv1_1,
		//		Dm_inv1_2, 0, 0,
		//		0, Dm_inv1_2, 0,
		//		0, 0, Dm_inv1_2;
		//	hessian += DF_DX2.transpose() * D2W_DFDF * DF_DX2;
		//}
		//else //Face vertex 3
		//{
		//	assert(iV == facePosVId(faceId, 2));
		//	force += F12.col(1);
		//	Eigen::Matrix < FloatingType, 6, 3> DF_DX3;
		//	DF_DX3 << Dm_inv2_1, 0, 0,
		//		0, Dm_inv2_1, 0,
		//		0, 0, Dm_inv2_1,
		//		Dm_inv2_2, 0, 0,
		//		0, Dm_inv2_2, 0,
		//		0, 0, Dm_inv2_2;
		//	hessian += DF_DX3.transpose() * D2W_DFDF * DF_DX3;
		//}

		FloatingType m1, m2;

		switch (faceVertexOrder)
		{
		case 0:
			m1 = -Dm_inv1_1 - Dm_inv2_1;
			m2 = -Dm_inv1_2 - Dm_inv2_2;
			force += -F12.col(0) - F12.col(1);
			break;
		case 1:
			m1 = Dm_inv1_1;
			m2 = Dm_inv1_2;
			force += F12.col(0);
			break;
		case 2:
			m1 = Dm_inv2_1;
			m2 = Dm_inv2_2;
			force += F12.col(1);
			break;
		default:
			break;
		}

		Mat3 h;
		assembleMembraneHessian(D2W_DFDF, m1, m2, h);
		hessian += h;
	}
	Mat3 dampingH = (hessian - old_hessian) * (objectParams().dampingStVK / pPhysicsParams->dt);
	Vec3 displacement = vertex(iV) - vertexPrevPos(iV);
	Vec3 dampingForce = dampingH * displacement;
	force -= dampingForce;
	hessian += dampingH;
}

void GAIA::VBDTriMeshStVk::accumlateBendingForceAndHessian(int vId, Vec3& force, Mat3& hessian)
{
	CFloatingType bendingStiffness = objectParams().bendingStiffness;
	CFloatingType dampingBending = objectParams().dampingBending;
	CFloatingType dt = pPhysicsParams->dt;
	size_t numberRelevantBendings = numRelevantBendings(vId);

	for (size_t iB = 0; iB < numberRelevantBendings; iB++)
	{
		const IdType eId = getVertexIthRelevantBending(vId, iB);
		const EdgeInfo& eInfo = getEdgeInfo(eId);
		const IdType edgeVertexOrder = getVertexIthRelevantBendingOrder(vId, iB);
		const Mat4& Q = edgeLaplacianQuadraticForms[eId];

		if (eInfo.fId2 == -1)
			// boundary edge
		{
			continue;
		}

		Eigen::Matrix<FloatingType, 4, 3> Xs;
		Xs.row(0) = vertex(eInfo.eV1);
		Xs.row(1) = vertex(eInfo.eV2);
		Xs.row(2) = vertex(eInfo.eV12Next);
		Xs.row(3) = vertex(eInfo.eV21Next);

		const Vec3 n1 = (Xs.row(1) - Xs.row(0)).cross(Xs.row(2) - Xs.row(0));
		const Vec3 n2 = (Xs.row(2) - Xs.row(3)).cross(Xs.row(1) - Xs.row(3));
		FloatingType n1_norm2 = n1.squaredNorm(), n2_norm2 = n2.squaredNorm();

		if (n1_norm2 < SQR(pPhysicsParams->degenerateTriangleThres)
			|| n2_norm2 < SQR(pPhysicsParams->degenerateTriangleThres))
			// degenerate triangle
		{
			continue;
		}

		const Eigen::Matrix<FloatingType, 4, 3> dE_dXs = bendingStiffness * Q * Xs;
		Mat3 hTemp = bendingStiffness * Q(edgeVertexOrder, edgeVertexOrder) * Mat3::Identity();
		force -= dE_dXs.row(edgeVertexOrder).transpose() + hTemp * (vertex(vId) - vertexPrevPos(vId)) * (dampingBending / dt);
		hTemp *= (1.f + dampingBending / dt);
		hessian += hTemp;

		//const Eigen::Matrix<FloatingType, 4, 3> dE_dXs = bendingStiffness * Q * Xs;
		//if (vId == eInfo.eV1)
		//{
		//	force += -dE_dXs.row(0);
		//	hessian += bendingStiffness * Q(0, 0) * Mat3::Identity();;
		//}
		//else if (vId == eInfo.eV2)
		//{
		//	force += -dE_dXs.row(1);
		//	hessian += bendingStiffness * Q(1, 1) * Mat3::Identity();;
		//}
		//else if (vId == eInfo.eV12Next) {
		//	force += -dE_dXs.row(2);
		//	hessian += bendingStiffness * Q(2, 2) * Mat3::Identity();
		//}
		//else if (vId == eInfo.eV21Next) {
		//	force += -dE_dXs.row(3);
		//	hessian += bendingStiffness * Q(3, 3) * Mat3::Identity();
		//}
		//else
		//{
		//	std::cerr << "Error! Inrelavant vertex encountered in the evaluation of the bending energy." << std::endl;
		//	getchar();
		//	exit(-1);
		//}
	}
}

FloatingType GAIA::VBDTriMeshStVk::backTracingLineSearchVBD(IdType vId, const Vec3& dx, FloatingType E0, FloatingType alpha, FloatingType c, FloatingType tau,
	int maxNumIters)
{
	FloatingType m = dx.squaredNorm();

	const Vec3 orgPos = vertex(vId);


	FloatingType bestAlpha = 0.f;
	FloatingType bestEnergy = E0;

	FloatingType meInertia = 0;
	FloatingType meElastic_stvk = 0;
	FloatingType meElastic_bending = 0;

	for (size_t iLineSearchIter = 0; iLineSearchIter < maxNumIters; iLineSearchIter++)
	{
		vertex(vId) = alpha * dx + orgPos;

		FloatingType e = evaluateVertexMeritEnergy(vId, meInertia, meElastic_stvk, meElastic_bending);

		if (e < bestEnergy)
		{
			bestAlpha = alpha;
			bestEnergy = e;
		}

		// first Wolfe condition 
		if (e < E0 - alpha * c * m)
		{
			break;
		}
		else
		{
			alpha = alpha * tau;
		}
	}

	return bestEnergy;
}


FloatingType GAIA::VBDTriMeshStVk::evaluatedAvgGradNorm()
{
	FloatingType gradNorm = 0.f;
	for (size_t iVert = 0; iVert < numVertices(); iVert++)
	{
		gradNorm += gradient.col(iVert).norm();
	}
	return gradNorm / numVertices();
}

FloatingType GAIA::VBDTriMeshStVk::evaluatedRegistrationEnergy()
{
	FloatingType registrationEnergy = 0.f;
	for (size_t iVert = 0; iVert < numVertices(); iVert++)
	{
		const TriMeshClosestPointQueryResult& targetPointResult = targetPointsQueyResults[iVert];
		if (targetPointResult.found)
		{
			registrationEnergy += 0.5f * pPhysicsParams->registrationEnergyWeight
				* (vertex(iVert) - targetPointResult.contactPts.back().closestPt).squaredNorm();
		}
	}
	return registrationEnergy;
}

FloatingType GAIA::VBDTriMeshStVk::evaluatedMeritEnergyElastic(FloatingType& meElastic_stvk, FloatingType& meElastic_bending)
{
	cpu_parallel_for(0, numFaces(), [&](int iFace) {
		Mat3x2 F;
		calculateDeformationGradient(iFace, F);
		Mat2 G = greenStrain(F);

		FloatingType e_uu = G(0, 0), e_vv = G(1, 1), e_uv = G(1, 0);

		meritEnergyElasticPerFace(iFace) = faceRestposeArea(iFace) * (0.5f * objectParams().lambda * pow(e_uu + e_vv, 2)
			+ objectParams().miu * (e_uu * e_uu + e_vv * e_vv + 2 * e_uv * e_uv));

		});

	// evaluateFaceNormals(true);
	cpu_parallel_for(0, numEdges(), [&](int iEdge) {
		const EdgeInfo& eInfo = pTopology->edgeInfos[iEdge];

		if (eInfo.fId2 == -1)
		{
			meritEnergyElasticPerEdge(iEdge) = 0.f;
			return;
		}

		Eigen::Matrix<FloatingType, 4, 3> Xs;
		Xs.row(0) = vertex(eInfo.eV1);
		Xs.row(1) = vertex(eInfo.eV2);
		Xs.row(2) = vertex(eInfo.eV12Next);
		Xs.row(3) = vertex(eInfo.eV21Next);

		const Vec3 n1 = (Xs.row(1) - Xs.row(0)).cross(Xs.row(2) - Xs.row(0));
		const Vec3 n2 = (Xs.row(2) - Xs.row(3)).cross(Xs.row(1) - Xs.row(3));
		FloatingType n1_norm = n1.norm(), n2_norm = n2.norm();

		if (n1_norm < pPhysicsParams->degenerateTriangleThres
			|| n2_norm < pPhysicsParams->degenerateTriangleThres)
		{
			meritEnergyElasticPerEdge(iEdge) = 0.f;
		}
		else
		{
			const Mat4& Q = edgeLaplacianQuadraticForms[iEdge];
			auto E = 0.5f * objectParams().bendingStiffness * (
				Xs.col(0).transpose() * Q * Xs.col(0)
				+ Xs.col(1).transpose() * Q * Xs.col(1)
				+ Xs.col(2).transpose() * Q * Xs.col(2));

			meritEnergyElasticPerEdge(iEdge) = E(0);
		}

		//Vec3 n1 = faceNormals.col(eInfo.fId1);
		//Vec3 n2 = faceNormals.col(eInfo.fId2);

		//const FloatingType cosTheta = n1.dot(n2);
		//meritEnergyElasticPerEdge(iEdge) = 0.25 * pObjectParamsMaterial->bendingStiffness * (1 - cosTheta);




		});

	meElastic_stvk = meritEnergyElasticPerFace.sum();
	meElastic_bending = meritEnergyElasticPerEdge.sum();

	FloatingType meritEnergyElastic = meElastic_stvk + meElastic_bending;
	if (!pPhysicsParams->associateGravityWithInertia)
	{
		for (size_t iV = 0; iV < numVertices(); iV++)
		{
			meritEnergyElastic += vertex(iV)(1) * pPhysicsParams->gravity(1) * vertexMass(iV);
		}
	}

	return meritEnergyElastic;
}

//void GAIA::VBDTriMeshStVk::accumulateMaterialGradient(const bool withInertia)
//{
//	evaluateInternalForce();
//
//	FloatingType dt = pPhysicsParams->dt;
//	FloatingType dtSqrReciprocal = 1.f / (dt * dt);
//	cpu_parallel_for(0, numVertices(), [&](int iV) {
//		// if not associated with inertia, gravity is associated as external force -(pPhysicsParams->gravity * vertexMass(iV))
//		Vec3 vgrad = -vertexExternalForces.col(iV) - vertexInternalForces.col(iV);
//
//		if (withInertia)
//		{
//			vgrad += vertexMass(iV) * (vertex(iV) - inertia.col(iV)) * (dtSqrReciprocal);
//		}
//
//		if (!pPhysicsParams->associateGravityWithInertia)
//		{
//			vgrad -= (pPhysicsParams->gravity * vertexMass(iV));
//		}
//
//		gradient.col(iV) += vgrad;
//		});
//
//}

int GAIA::VBDTriMeshStVk::tearMesh(IdType v1, IdType v2)
{
	int ret = GAIA::TriMeshFEM::tearMesh(v1, v2);
	vertexExternalForces.conservativeResize(Eigen::NoChange, numVertices());
	if (ret == 1) {
		vertexExternalForces.col(numVertices() - 1) = vertexExternalForces.col(v1);
	}
	else if (ret == 2) {
		vertexExternalForces.col(numVertices() - 1) = vertexExternalForces.col(v2);
	}
	else if (ret == 3) {
		vertexExternalForces.col(numVertices() - 2) = vertexExternalForces.col(v1);
		vertexExternalForces.col(numVertices() - 1) = vertexExternalForces.col(v2);
	}
	return ret;
}



bool GAIA::VBDObjectParamsTriMeshStVk::fromJson(nlohmann::json& objectParam)
{
	VBDObjectParamsTriMesh::fromJson(objectParam);

	EXTRACT_FROM_JSON(objectParam, lambda);
	EXTRACT_FROM_JSON(objectParam, miu);
	EXTRACT_FROM_JSON(objectParam, dampingStVK);
	EXTRACT_FROM_JSON(objectParam, dampingBending);
	return true;
}

bool GAIA::VBDObjectParamsTriMeshStVk::toJson(nlohmann::json& objectParam)
{
	VBDObjectParamsTriMesh::toJson(objectParam);

	PUT_TO_JSON(objectParam, lambda);
	PUT_TO_JSON(objectParam, miu);
	PUT_TO_JSON(objectParam, dampingStVK);
	PUT_TO_JSON(objectParam, dampingBending);
	return true;
}



#pragma once
#include <functional>
#include <memory>

#include "../Deformer.h"
#include "../PBDPhysics.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

namespace EBD {

	struct DeformerRigid1DoF : public BaseDeformer {
		DeformerRigid1DoF(const Vec3 inO, const Vec3 inAxis, double inAngularVelocity, PBDPhysics& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inRotationEndTime = -1) :
			o(inO), axis(inAxis), angularVelocity(inAngularVelocity), selectedMeshes(inSelectedMeshes),
			selectedVerts(inSelectedVerts), rotationEnd(inRotationEndTime)
		{
			axis = axis / axis.norm();
			allPts = 0;
			for (int iMesh = 0; iMesh <  selectedMeshes.size(); ++iMesh)
			{
				allPts += selectedVerts[iMesh].size();
			}
			orgPs.resize(POINT_VEC_DIMS, allPts);
			int ivAll = 0;

			for (int iMesh : selectedMeshes)
			{
				int meshId = selectedMeshes[iMesh];
				TetMeshFEM* pTM = physics.tMeshes[meshId].get();
				for (int iV = 0; iV < selectedVerts[iMesh].size(); ++iV) {
					int  vId = selectedVerts[iMesh][iV];
					pTM->verticesCollisionDetectionEnabled[vId] = false;

					orgPs.col(ivAll) = pTM->vertex(vId);
					++ivAll;
				}
			}
		}

		virtual void operator()(PBDPhysics* physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
			TVerticesMat newPos;
			int ivAll = 0;

			for (int iMesh = 0; iMesh < selectedMeshes.size(); ++iMesh)
			{
				int meshId = selectedMeshes[iMesh];
				TetMeshFEM* pTM = physics->tMeshes[meshId].get();
				newPos.resizeLike(orgPs);
				for (int iV = 0; iV < selectedVerts[iMesh].size(); ++iV) {
					int  vId = selectedVerts[iMesh][iV];
					// compute rigid deformation
					newPos.col(ivAll) = pTM->vertex(vId);
					//std::cout << "selected pt new pos: " << newPos.col(iV).transpose() << std::endl;
					++ivAll;
				}
			}

			Eigen::Matrix<double, Eigen::Dynamic, POINT_VEC_DIMS> A = orgPs.cast<double>().transpose(); 
			A = (A.rowwise() - o.cast<double>().transpose());
			Eigen::Matrix<double, Eigen::Dynamic, POINT_VEC_DIMS> B = newPos.cast<double>().transpose();
			B = (B.rowwise() - o.cast<double>().transpose());
			Eigen::Matrix<double, 1, 3> centroid_A = A.colwise().mean();
			Eigen::Matrix<double, 1, 3> centroid_B = B.colwise().mean();
			Eigen::Matrix<double, Eigen::Dynamic, POINT_VEC_DIMS> AA = (A.rowwise() - centroid_A);
			Eigen::Matrix<double, Eigen::Dynamic, POINT_VEC_DIMS> BB = (B.rowwise() - centroid_B);

			Mat3d H = AA.transpose() * BB;
			// Eigen::JacobiSVD<Mat3d, Eigen::ComputeFullU | Eigen::ComputeFullV> svd(H);
			Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Vec3d S = svd.singularValues();
			Mat3d U = svd.matrixU();
			Mat3d V = svd.matrixV();
			Mat3d R = V * U.transpose();
			// Mat3d R = V * U.transpose();

			if (R.determinant() < 0) {
				V.block<3, 1>(0, 2) *= -1;
				R = V * U.transpose();
			}

			Mat3d HRecover = U * S.asDiagonal() * V.transpose();

			Eigen::AngleAxisf angleAxis(R.cast<float>());

			FloatingType rotAngle = angleAxis.angle() * axis.dot(angleAxis.axis());

			angleAxis.angle() = rotAngle;
			angleAxis.axis() = axis;

			TVerticesMat AA_rot = angleAxis.toRotationMatrix() * (AA.cast<float>().transpose());

			ivAll = 0;
			for (int iMesh = 0; iMesh < selectedMeshes.size(); ++iMesh)
			{
				int meshId = selectedMeshes[iMesh];
				TetMeshFEM* pTM = physics->tMeshes[meshId].get();
				for (int iV = 0; iV < selectedVerts[iMesh].size(); ++iV) {
					int  vId = selectedVerts[iMesh][iV];
					// compute rigid deformation
					pTM->vertex(vId) = AA_rot.col(ivAll) + centroid_A.cast<float>().transpose() + o;
					++ivAll;
				}
			}

		};

		Vec3 o;
		Vec3 axis;

		double angularVelocity;

		double rotationEnd;

		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		TVerticesMat orgPs;
		int allPts = -1;
	};
}
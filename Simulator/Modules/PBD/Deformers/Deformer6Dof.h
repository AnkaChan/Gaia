#pragma once
#include <functional>
#include <memory>

#include "../Deformer.h"
#include "../PBDPhysics.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

namespace EBD {
	struct Deformer6Dof : public BaseDeformer {
		Deformer6Dof(const Vec3 inO, const Vec3 inAxis, const Vec3 targetTranslation , double inTargetRot, 
			PBDPhysics& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inRotationEndTime = -1) :
			rotCenter(inO),
			axis(inAxis), 
			targetRot(inTargetRot), 
			targetTranslation(targetTranslation),
			selectedMeshes(inSelectedMeshes),
			selectedVerts(inSelectedVerts), 
			deformEnd(inRotationEndTime)
		{
			for (int iMesh : selectedMeshes)
			{
				orgOPs.emplace_back();
				TetMeshFEM* pTM = physics.tMeshes[iMesh].get();
				for (int iV : selectedVerts[iMesh]) {
					Vec3 p = pTM->vertex(iV);
					pTM->verticesCollisionDetectionEnabled[iV] = false;

					Vec3 P = { p[0], p[1], p[2] };
					Vec3 OP = P - rotCenter;

					orgOPs.back().push_back(OP);

				}
			}
		}

		virtual void operator()(PBDPhysics* physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
			if (deformEnd <= 0)
			{
				return;
			}
			if (curTime >= deformEnd && deformEnd > 0)
			{
				curTime = deformEnd;
			}
			double w = curTime / deformEnd;
			double theta = w * targetRot;
		

			Eigen::AngleAxisf angleAxis;

			angleAxis.angle() = theta;
			angleAxis.axis() = axis;

			Eigen::Matrix3f R;
			R = angleAxis.toRotationMatrix();
			// Eigen::Vector3d O = { o[0], o[1], o[2] };
			// Eigen::Vector3d N = { axis[0], axis[1], axis[2] };

			for (int iMesh : selectedMeshes)
			{
				TetMeshFEM* pTM = physics->tMeshes[iMesh].get();
				int iV = 0;
				for (int vId : selectedVerts[iMesh]) {
					Vec3Block p = pTM->vertex(vId);

					const Vec3& OP = orgOPs[iMesh][iV];

					const Vec3 RootPRot = R * OP;
					Vec3 PRot = rotCenter + RootPRot;

					// Eigen::Vector3d PRot = root + RootP;
					//std::cout << "Rotate from: " << p.transpose() << " to " << PRot.transpose() << std::endl;

					PRot = PRot + w * targetTranslation;
					
					p[0] = PRot[0];
					p[1] = PRot[1];
					p[2] = PRot[2];
					++iV;
				}
			}

		};

		Vec3 targetTranslation;
		Vec3 axis;
		Vec3 rotCenter;
		double targetRot;
		double deformEnd;
		
		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		std::vector<std::vector<Vec3>> orgOPs;
	public:
		EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	};
}
#pragma once
#include <functional>
#include <memory>

#include "../Deformer.h"
#include "../PBDPhysics.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

namespace EBD {
	struct DeformerRotator : public BaseDeformer {
		DeformerRotator(const Vec3 inO, const Vec3 inAxis, double inAngularVelocity, PBDPhysics& physics, const std::vector<int>& inSelectedMeshes,
			const std::vector<std::vector<int>>& inSelectedVerts, double inRotationEndTime = -1) :
			o(inO), axis(inAxis), angularVelocity(inAngularVelocity), selectedMeshes(inSelectedMeshes), 
			selectedVerts(inSelectedVerts), rotationEnd(inRotationEndTime)
		{
			for (int iMesh : selectedMeshes)
			{
				TetMeshFEM* pTM = physics.tMeshes[iMesh].get();
				roots.emplace_back();
				RootPs.emplace_back();
				for (int iV : selectedVerts[iMesh]) {
					Vec3 p = pTM->vertex(iV);
					pTM->verticesCollisionDetectionEnabled[iV] = false;

					Vec3 P = { p[0], p[1], p[2] };
					Vec3 OP = P - o;
					Vec3 root = OP.dot(axis) * axis;
					Vec3 RootP = P - root;

					roots.back().push_back(root);
					RootPs.back().push_back(RootP);

				}
			}
		}

		virtual void operator()(PBDPhysics* physics, double curTime, int iFrame, int iSubstep, int iIter, double dt) {

			if (curTime >= rotationEnd && rotationEnd > 0)
			{
				curTime = rotationEnd;
			}
			double theta = curTime * angularVelocity;


			double ux = axis[0];
			double uy = axis[1];
			double uz = axis[2];

			Eigen::Matrix3f R;
			R << cos(theta) + ux * ux * (1 - cos(theta)), ux* uy* (1 - cos(theta)) - uz * sin(theta), ux* uz* (1 - cos(theta)) + uy * sin(theta),
				uy* ux* (1 - cos(theta)) + uz * sin(theta), cos(theta) + uy * uy * (1 - cos(theta)), uy* uz* (1 - cos(theta)) - ux * sin(theta),
				uz* ux* (1 - cos(theta)) - uy * sin(theta), uz* uy* (1 - cos(theta)) + ux * sin(theta), cos(theta) + uz * uz * (1 - cos(theta));

			// Eigen::Vector3d O = { o[0], o[1], o[2] };
			// Eigen::Vector3d N = { axis[0], axis[1], axis[2] };

			for (int iMesh : selectedMeshes)
			{
				TetMeshFEM* pTM = physics->tMeshes[iMesh].get();
				int iV = 0;
				for (int vId : selectedVerts[iMesh]) {
					Vec3Block p = pTM->vertex(vId);

					const Vec3& root = roots[iMesh][iV];

					const Vec3& RootP = RootPs[iMesh][iV];

					const Vec3 RootPRot = R * RootP;
					const Vec3 PRot = root + RootPRot;

					// Eigen::Vector3d PRot = root + RootP;
					//std::cout << "Rotate from: " << p.transpose() << " to " << PRot.transpose() << std::endl;
					
					p[0] = PRot[0];
					p[1] = PRot[1];
					p[2] = PRot[2];
					++iV;
				}
			}

		};

		Vec3 o;
		Vec3 axis;

		double angularVelocity;

		double rotationEnd;

		std::vector<int> selectedMeshes;
		std::vector<std::vector<int>> selectedVerts;
		std::vector<std::vector<Vec3>> roots;
		std::vector<std::vector<Vec3>> RootPs;
	};
}
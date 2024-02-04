#pragma once
#include <functional>
#include <memory>

#include "Deformer.h"
#include "../SoftBody/SoftBody.h"

#include <cmath>
#include <Eigen/StdVector>
#include "Eigen/core"

struct DeformerRigidPuller : public BaseDeformer {
	DeformerRigidPuller(CPoint inVelocity, SoftBodyManager* pSoftBodyManager, const std::vector<int> & inSelectedMeshes, const std::vector<std::vector<int>>& inSelectedVerts, double in_endTime=-1) :
		velocity(inVelocity), selectedMeshes(inSelectedMeshes), selectedVerts(inSelectedVerts), endTime(in_endTime)
	{

		for (int iMesh : selectedMeshes)
		{
			orgPos.emplace_back();
			for (int iV : selectedVerts[iMesh]) {
				CPoint& p = (*pSoftBodyManager->softBodies[iMesh]->mesh.verts[iV]);

				orgPos.back().push_back(p);
			}
		}
	}

	virtual void operator()(SoftBodyManager* pSoftBodyManager, double curTime, int iFrame, int iSubstep, int iIter, double dt) {

		for (int iMesh : selectedMeshes)
		{
			int iV = 0;
			for (int vId : selectedVerts[iMesh]) {
				CPoint& p = (*pSoftBodyManager->softBodies[iMesh]->mesh.verts[vId]);

				if (endTime > 0 && curTime > endTime)
				{
					curTime = endTime;
				}
				CPoint newP = orgPos[iMesh][iV] + velocity * curTime;

				p = newP;
				++iV;
			}
		}
	};

	CPoint velocity;
	
	std::vector<int> selectedMeshes;
	std::vector<std::vector<int>> selectedVerts;
	std::vector<std::vector<CPoint>> orgPos;

	double endTime = -1;
};
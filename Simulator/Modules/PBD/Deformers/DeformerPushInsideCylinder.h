#pragma once
#include <functional>
#include <memory>

#include "Deformer.h"
#include "../SoftBody/SoftBody.h"

struct DeformerPushToCylinder : public BaseDeformer {
	DeformerPushToCylinder(CPoint inCylinerCenter, double inRadius, double inAddVelocity) :
		cylinerCenter(inCylinerCenter), radius(inRadius), addVelocity(inAddVelocity)
	{

	}

	virtual void operator()(SoftBodyManager* pSoftBodyManager, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
		if (iIter == 0)
		{
			for (size_t iMesh = 0; iMesh < pSoftBodyManager->softBodies.size(); iMesh++)
			{
				SoftBody& softbody = *pSoftBodyManager->softBodies[iMesh];
				for (size_t iV = 0; iV < softbody.mesh.verts.size() ; iV++)
				{
					CPoint& p = *softbody.mesh.verts[iV];
					CPoint diff = cylinerCenter - p;
					// assume the center axis of the cylinder is (0, 1, 0)
					diff[1] = 0;
					double disToCenter = diff.norm();

					CPoint diffN = diff / diff.norm();
					if (disToCenter > radius)
					{
						p += diffN * addVelocity * dt;
					}
				}
			}
			
		}
	};

	CPoint cylinerCenter;
	double radius;
	double addVelocity;
};
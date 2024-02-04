#include "Deformer.h"

struct PostDeformerStablizer : public BaseDeformer {
	PostDeformerStablizer(int inDampingStartFrame, double inPerSecondDamping=0.997) :
		dampingStartFrame(inDampingStartFrame), perSecondDamping(inPerSecondDamping)
	{
	}

	virtual void operator()(SoftBodyManager* pSoftBodyManager, double curTime, int iFrame, int iSubstep, int iIter, double dt) {
		if (iFrame >= dampingStartFrame) {
			
			double stepInvariantVelDamping = std::pow(perSecondDamping, (pSoftBodyManager->curPhysics.timeStep / pSoftBodyManager->curPhysics.numSubsteps));
			
			int framesSinceDampingStart = iFrame - dampingStartFrame;
			double dampingRate = pow(stepInvariantVelDamping, framesSinceDampingStart * pSoftBodyManager->curPhysics.timeStep);
			std::cout << dampingRate << "\n";

			for (size_t iSb = 0; iSb < pSoftBodyManager->softBodies.size(); iSb++)
			{
				SoftBody& sb = *pSoftBodyManager->softBodies[iSb];
				for (size_t iVert = 0; iVert < sb.mesh.verts.size(); iVert++)
				{
					CPoint& prevPos = sb.mesh.prevPos[iVert];
					*sb.mesh.verts[iVert] = prevPos + (*sb.mesh.verts[iVert] - prevPos) * dampingRate;
				}
			}
		}
	}


	int dampingStartFrame;
	double perSecondDamping = 0.997;

};
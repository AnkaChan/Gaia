#include "../Framework/BasePhysicsFramework.h"
#include "VBD_BaseMaterial.h"
#include "../Parallelization/CPUParallelization.h"

using namespace GAIA;

void VBDBaseTetMesh::forwardStep()
{
	CFloatingType dt = pPhysicsParams->dt;
	CFloatingType gravNorm = pPhysicsParams->gravity.norm();
	Vec3 gravDir = pPhysicsParams->gravity / gravNorm;
#ifdef OUTPUT_INITIALIZATION_GRAV_NORM
	std::vector<FloatingType> gravNormsForInitialization;
#endif // OUTPUT_INITIALIZATION_GRAV_NORM

	switch (pObjParamsVBD->initializationType)
	{
	case 0:
		// full inertia
		positions() = inertia;
		break;
	case 1:
		// previous posistion
		// positions() = positionsPrev();
		break;

	case 2:

		positions() = positionsPrev() + dt * velocitiesPrev();
		break;

	case 3:
		positions() = 0.5f * positions() + 0.5f * inertia;

		break;

	case 4:
		positions() = positions() + dt * velocitiesPrev();
		positions().colwise() += dt * dt * pObjParamsVBD->initRatio_g * pPhysicsParams->gravity;

		break;
	case 5:
		forwardStepSimplecticEuler();
		break;

	case 6:
		// adaptive 

		for (int iV = 0; iV < numVertices(); iV++)
		{
			FloatingType accelerationComponent = 0.f;
			if (hasApproxAcceleration)
			{
				accelerationComponent = acceletration.col(iV).dot(gravDir);
				accelerationComponent = accelerationComponent < gravNorm ? accelerationComponent : gravNorm;
				accelerationComponent = accelerationComponent > 1e-5f ? accelerationComponent : 0.f;
				mVertPos.col(iV) = positionsPrev().col(iV) + dt * velocitiesPrev().col(iV) + dt * dt * gravDir * accelerationComponent;
			}
			else
			{
				accelerationComponent = gravNorm;
				//positions() = inertia;
			}
#ifdef OUTPUT_INITIALIZATION_GRAV_NORM
			gravNormsForInitialization.push_back(accelerationComponent);
#endif // OUTPUT_INITIALIZATION_GRAV_NORM
		}

#ifdef OUTPUT_INITIALIZATION_GRAV_NORM
		{
			if (pPhysicsFramework->substep==0)
			{
				std::vector<FloatingType> gravNormsForInitializationSurfaceVerts;

				for (size_t iSurfaceV = 0; iSurfaceV < surfaceVIds().size(); iSurfaceV++)
				{
					gravNormsForInitializationSurfaceVerts.push_back(gravNormsForInitialization[surfaceVIds()(iSurfaceV)]);
				}

				std::string debugFolder = pPhysicsFramework->getDebugFolder();

				std::ostringstream aSs;
				aSs << debugFolder << "/InitializationGravNorm_Mesh_"
					<< std::setfill('0') << std::setw(6) << meshId
					<< "_Frame_" << std::setfill('0') << std::setw(8) << pPhysicsFramework->frameId
					//<< "Substep_" << std::setfill('0') << std::setw(4) << pPhysicsFramework->substep
					<< ".json";

				nlohmann::json j;
				j["InitializationGravNorm"] = gravNormsForInitializationSurfaceVerts;
				std::string filePath = aSs.str();
				bool retVal = MF::saveJson(filePath, j);
			}
		}
#endif // OUTPUT_INITIALIZATION_GRAV_NORM
		/*for (int iV = 0; iV < numVertices(); iV++)
		{
			if (hasApproxAcceleration)
			{

				FloatingType accelerationNorm = acceletration.col(iV).norm();
				if (accelerationNorm > 1e-5f)
				{
					Vec3 accelerationDir = acceletration.col(iV) / accelerationNorm;
					FloatingType gravComponent = accelerationDir.dot(pPhysicsParams->gravity);

					if (gravComponent < 0)
					{
						gravComponent = 0.f;
					}else if(gravComponent > accelerationNorm) {
						gravComponent = accelerationNorm;
					}
					mVertPos.col(iV) = positionsPrev().col(iV) + dt * velocitiesPrev().col(iV) + dt * dt * accelerationDir * gravComponent;
				}
				else
				{
					mVertPos.col(iV) = positionsPrev().col(iV) + dt * velocitiesPrev().col(iV) ;
				}
			}
			else
			{
				positions() = inertia;
			}

		}*/
		//mVertPos = positionsPrev() + dt * velocitiesPrev() + dt * dt * acceletration;

		break;
	default:
		break;
	}
}

void VBDBaseTetMesh::forwardStepSimplecticEuler()
{
	CFloatingType dt = pPhysicsParams->dt;
	evaluateInternalForce();
	evaluateExternalForce();

	cpu_parallel_for(0, numVertices(), [&](int iV)
		{
			mVelocity.col(iV) = mVelocitiesPrev.col(iV) + dt * (vertexExternalForces.col(iV) + vertexInternalForces.col(iV))
				* vertexInvMass(iV);

			mVertPos.col(iV) += mVelocity.col(iV) * dt;
		});
}
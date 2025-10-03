#include "VBDBaseTriMesh.h"
#include <Parallelization/CPUParallelization.h>

using namespace GAIA;

bool GAIA::VBDObjectParamsTriMesh::fromJson(nlohmann::json& objectParam)
{
	TriMeshParams::fromJson(objectParam);

	EXTRACT_FROM_JSON(objectParam, exponentialVelDamping);
	EXTRACT_FROM_JSON(objectParam, constantVelDamping);
	EXTRACT_FROM_JSON(objectParam, initializationType);

	return true;
}

bool GAIA::VBDObjectParamsTriMesh::toJson(nlohmann::json& objectParam)
{
	TriMeshParams::toJson(objectParam);

	PUT_TO_JSON(objectParam, exponentialVelDamping);
	PUT_TO_JSON(objectParam, constantVelDamping);
	PUT_TO_JSON(objectParam, initializationType);

	return true;
}

void GAIA::VBDBaseTriMesh::initialize(TriMeshParams::SharedPtr inObjectParams, BaseClothPhsicsFramework::Ptr in_pPhysics)
{
	TriMeshFEM::initialize(inObjectParams, true);
	pPhysics = in_pPhysics;
	pPhysicsParams = (VBDClothPhysicsParameters::Ptr)pPhysics->basePhysicsParams.get();

	positionsNext.resize(3, numVertices());

	edgeDegenrateMask.resize(numEdges());

	triangleInternalForces.resize(9, numFaces());
	edgeInternalForces.resize(12, pTopology->numEdges);

	vertexInternalForces.resize(3, numVertices());
	vertexExternalForces.resize(3, numVertices());
	inertia.resize(3, numVertices());
	meritEnergy_Inertia_PerVertex.resize(numVertices());
	meritEnergyElasticPerFace.resize(numFaces());
	meritEnergyElasticPerEdge.resize(pTopology->numEdges);
	gradient.resize(3, numVertices());
	gradient.setZero();

	positionsAtPrevCollisionDetection.resize(3, numVertices());
	positionsAtPrevCollisionDetection.setZero();
	vertexConvervativeBounds.resize(numVertices());

	vertexConstraints.resize(numVertices());
	activeCollisionMask.resize(numVertices());
}

void GAIA::VBDBaseTriMesh::forwardStep_intertia()
{
	FloatingType dt = pPhysicsParams->dt;
	evaluateExternalForce();
	for (size_t iV = 0; iV < numVertices(); iV++)
	{
		velocities.col(iV) += dt * vertexExternalForces.col(iV) * vertexInvMass(iV);
	}
	if (pPhysicsParams->associateGravityWithInertia)
	{
		velocities.colwise() += dt * pPhysicsParams->gravity;
	}
	handleVeclocityConstraint();
	inertia = positions() + dt * velocities;

	CFloatingType gravNorm = pPhysicsParams->gravity.norm();
	Vec3 gravDir = pPhysicsParams->gravity / gravNorm;
	switch (objectParamsBase().initializationType)
	{
	case 0:
		// full inertia
		positions() = inertia;
		break;
	case 1:
		// adaptive 

		for (int iV = 0; iV < numVertices(); iV++)
		{
			FloatingType accelerationComponent = 0.f;
			if (hasApproxAcceleration)
			{
				accelerationComponent = acceletration.col(iV).dot(gravDir);
				accelerationComponent = accelerationComponent < gravNorm ? accelerationComponent : gravNorm;
				accelerationComponent = accelerationComponent > 1e-5f ? accelerationComponent : 0.f;
				positions().col(iV) = positionsPrev.col(iV) + dt * velocitiesPrev.col(iV) + dt * dt * gravDir * accelerationComponent;
			}
			else
			{
				accelerationComponent = gravNorm;
				positions().col(iV) = inertia.col(iV);
			}
#ifdef OUTPUT_INITIALIZATION_GRAV_NORM
			gravNormsForInitialization.push_back(accelerationComponent);
#endif // OUTPUT_INITIALIZATION_GRAV_NORM
		}

#ifdef OUTPUT_INITIALIZATION_GRAV_NORM
		{
			if (pPhysicsFramework->substep == 0)
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
		break;

	case 2:
		positions() = positionsPrev;
		break;

	case 3:
		// inertia within conservative bounds
		for (int iV = 0; iV < numVertices(); iV++)
		{
			CFloatingType convservativeBounds = vertexConvervativeBounds[iV];
			Vec3 inertiaStep = (inertia.col(iV) - positionsPrev.col(iV));
			FloatingType inertiaSize = inertiaStep.norm();

			if (pPhysicsParams->handleCollision && inertiaSize > convservativeBounds && inertiaSize)
				// moved out of the conservative bounds, truncation needed, and collision detection is required
			{
				FloatingType accumulatedDisplacementSizeNew = convservativeBounds;
				inertiaStep = inertiaStep * (accumulatedDisplacementSizeNew / inertiaSize);
			}
			vertex(iV) = positionsPrev.col(iV) + inertiaStep;
		}

		break;

	case 4:
		break;
	default:
		break;
	}
}

void GAIA::VBDBaseTriMesh::applyInitialStep()
{
	positionsPrev = positions();
	if (hasVelocitiesPrev)
	{
		acceletration = (velocities - velocitiesPrev) / pPhysicsParams->dt;
		hasApproxAcceleration = true;
	}
	else
	{
		hasVelocitiesPrev = true;
	}
	velocitiesPrev = velocities;
	forwardStep_intertia();
	evaluateExternalForce();
}

//void GAIA::VBDBaseTriMesh::forwardStep_symplecticEuler()
//{
//	positionsPrev = positions();
//	velocitiesPrev = velocities;
//
//	evaluateInternalForce();
//	evaluateExternalForce();
//	handleInternalForceConstraint();
//	FloatingType dt = pPhysicsParams->dt;
//	for (size_t iV = 0; iV < numVertices(); iV++)
//	{
//		velocities.col(iV) += dt * (vertexExternalForces.col(iV) + vertexInternalForces.col(iV)) * vertexInvMass(iV);
//	}
//
//	// std::cout << "vertexInternalForces:\n" << vertexInternalForces.transpose() << std::endl;
//
//	velocities.colwise() += dt * pPhysicsParams->gravity;
//	handleVeclocityConstraint();
//
//	positions() += dt * velocities;
//	// std::cout << "positions:\n" << positions() << std::endl;
//
//
//}

void GAIA::VBDBaseTriMesh::evaluateExternalForce()
{
	vertexExternalForces.setZero();
}

void GAIA::VBDBaseTriMesh::clearGradient()
{
	gradient.setZero();
	vertexInternalForces.setZero();
}

FloatingType GAIA::VBDBaseTriMesh::evaluatedMeritEnergy(FloatingType& meInertia, FloatingType& meElastic_stvk, FloatingType& meElastic_bending)
{
	meInertia = evaluatedMeritEnergyInertia();
	FloatingType meElastic = evaluatedMeritEnergyElastic(meElastic_stvk, meElastic_bending);
	return meElastic + meInertia;
}

FloatingType GAIA::VBDBaseTriMesh::evaluatedMeritEnergyInertia()
{
	FloatingType dt = pPhysicsParams->dt;
	FloatingType dtSqrReciprocal = 1.f / (dt * dt);
	cpu_parallel_for(0, numVertices(), [&](int iV) {
		meritEnergy_Inertia_PerVertex(iV) = 0.5f * vertexMass(iV) * (vertex(iV) - inertia.col(iV)).squaredNorm() * (dtSqrReciprocal);
		});

	return meritEnergy_Inertia_PerVertex.sum();
}

//void GAIA::VBDBaseTriMesh::forwardStep_IE_GD()
//{
//	positionsPrev = positions();
//	velocitiesPrev = velocities;
//
//	forwardStep_intertia();
//	evaluateExternalForce();
//
//	FloatingType dt = pPhysicsParams->dt;
//	FloatingType stepSize = pPhysicsParams->stepSize;
//	FloatingType dtSqrReciprocal = 1.f / (dt * dt);
//
//	FloatingType meInertia, meElastic_stvk, meElastic_bending;
//	FloatingType me = evaluatedMeritEnergy(meInertia, meElastic_stvk, meElastic_bending);
//
//	std::cout << "Initial merit energy: " << me << "\n";
//	for (size_t iteration = 0; iteration < pPhysicsParams->maxNumIteration; iteration++)
//	{
//		evaluateInternalForce();
//		handleInternalForceConstraint();
//
//		if (iteration == 0)
//		{
//			// backtracing line search to stop 
//		}
//
//		cpu_parallel_for(0, numVertices(), [&](int iV) {
//			Vec3 VInertia = inertia.col(iV);
//			Vec3 vgrad = vertexMass(iV) * (vertex(iV) - VInertia) * (dtSqrReciprocal)
//				-vertexExternalForces.col(iV) - vertexInternalForces.col(iV); // gravity is associated as external force -(pPhysicsParams->gravity * vertexMass(iV))
//
//			if (!pPhysicsParams->associateGravityWithInertia)
//			{
//				vgrad -= (pPhysicsParams->gravity * vertexMass(iV));
//			}
//
//				vertex(iV) -= vgrad * stepSize;
//			});
//
//		if (!(iteration % pPhysicsParams->evaluationSteps))
//		{
//			me = evaluatedMeritEnergy(meInertia, meElastic_stvk, meElastic_bending);
//			std::cout << "Merit energy at step " << pPhysics->substep << " iter " << iteration << ": "
//				<< me << " | meInertia: " << meInertia << " | energy_stvk: " << meElastic_stvk << " | energy_bending: " << meElastic_bending << "\n";
//		}
//
//	}
//	// !!! dont forget to update velocity!!
//	velocities = (positions() - positionsPrev) / dt;
//}

//void GAIA::VBDBaseTriMesh::forwardStep_IE_VBD()
//{
//	positionsPrev = positions();
//	velocitiesPrev = velocities;
//
//	forwardStep_intertia();
//	evaluateExternalForce();
//
//	FloatingType dt = pPhysicsParams->dt;
//
//	FloatingType meInertia, meElastic_stvk, meElastic_bending;
//	FloatingType me = evaluatedMeritEnergy(meInertia, meElastic_stvk, meElastic_bending);
//
//	std::cout << "Initial merit energy: " << me << "\n";
//	for (size_t iteration = 0; iteration < pPhysicsParams->maxNumIteration; iteration++)
//	{
//		for (size_t iV = 0; iV < numVertices(); iV++)
//		{
//			VBDStep(iV);
//		}
//
//	}
//	// !!! dont forget to update velocity!!
//	velocities = (positions() - positionsPrev) / dt;
//}

void GAIA::VBDBaseTriMesh::accumlateInertiaForceAndHessian(int iV, Vec3& force, Mat3& hessian)
{
	CFloatingType dtSqrReciprocal = pPhysicsParams->dtSqrReciprocal;

	force += vertexMass(iV) * (inertia.col(iV) - vertex(iV)) * (dtSqrReciprocal);

	hessian += Mat3::Identity() * vertexMass(iV) * (dtSqrReciprocal);
}

void GAIA::VBDBaseTriMesh::accumlateSewingForceAndHessian(int iV, Vec3& force, Mat3& hessian)
{
	const auto& sewingStiffness = pObjectParams->sewingStiffness;
	for (const auto& p : sewingVMap[iV]) {
		int sid = p.first;
		int order = p.second;
		FloatingType m = 1.f;
		FloatingType lbda = sewingRatio[sid];
		if (order == 1) {
			m = -(1 - lbda);
		}
		else if (order == 2) {
			m = -lbda;
		}
		int vi0 = sewingVertices[sid][0];
		int vi1 = sewingVertices[sid][1];
		int vi2 = sewingVertices[sid][2];
		Vec3 x0 = vertex(vi0);
		Vec3 x1 = vertex(vi1);
		Vec3 x2 = vertex(vi2);
		Vec3 x = x0 - (1 - lbda) * x1 - lbda * x2;
		force -= (m * sewingStiffness) * x;
		hessian += (m * m * sewingStiffness) * Mat3::Identity();
	}
}

void GAIA::VBDBaseTriMesh::handleInternalForceConstraint()
{
	for (size_t iFixedV = 0; iFixedV < pObjectParams->fixedPoints.size(); iFixedV++)
	{
		vertexInternalForces.col(pObjectParams->fixedPoints[iFixedV]).setZero();
	}
}

void GAIA::VBDBaseTriMesh::handleVeclocityConstraint()
{
	for (size_t iFixedV = 0; iFixedV < pObjectParams->fixedPoints.size(); iFixedV++)
	{
		velocities.col(pObjectParams->fixedPoints[iFixedV]).setZero();
	}
}
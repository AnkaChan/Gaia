#include "PBDTetMeshFEM.h"
#include "PBDPhysics.h"

bool GAIA::ObjectParamsPBD::fromJson(nlohmann::json& objectParam)
{
	ObjectParams::fromJson(objectParam);

	EXTRACT_FROM_JSON(objectParam, exponentialVelDamping);
	EXTRACT_FROM_JSON(objectParam, constantVelDamping);
	EXTRACT_FROM_JSON(objectParam, friction_dynamic);
	EXTRACT_FROM_JSON(objectParam, friction_static);

	return true;
}

bool GAIA::ObjectParamsPBD::toJson(nlohmann::json& objectParam)
{
	ObjectParams::toJson(objectParam);

	PUT_TO_JSON(objectParam, exponentialVelDamping);
	PUT_TO_JSON(objectParam, constantVelDamping);
	PUT_TO_JSON(objectParam, friction_dynamic);
	PUT_TO_JSON(objectParam, friction_static);

	return true;
}


void GAIA::PBDTetMeshFEM::initialize(ObjectParams::SharedPtr inMaterialParams, TetMeshMF::SharedPtr pTM_MF, PBDPhysics* inPPBDPhysics)
{
	TetMeshFEM::initialize(inMaterialParams, pTM_MF);

	std::cout << "Added tetmesh: " << inMaterialParams->path << "\n"
		<< "with " << numVertices() << " vertices and " << numTets() << "tets.\n";
	pPBDPhysics = inPPBDPhysics;
	dt = inPPBDPhysics->dt;
}

void GAIA::PBDTetMeshFEM::applyInitialGuess()
{
	// vel.vel[i] = (mesh.vel[i]) + ((gravity)*dt);
	mVertPrevPos = mVertPos;

	if (pObjectParams->hasNoGravZone)
	{
		for (int iVert = 0; iVert < m_nVertices; iVert++) {
			if (mVertPos(GRAVITY_AXIS, iVert) <= pObjectParams->noGravZoneThreshold)
			{
				mVelocity.col(iVert) += pPBDPhysics->physicsParams().gravity * pPBDPhysics->dt;
			}
		}
	}
	else {
		mVelocity.colwise() += pPBDPhysics->physicsParams().gravity * pPBDPhysics->dt;
	}
	
	mVertPos += mVelocity * pPBDPhysics->dt;
	inertia = mVertPos;
}

void GAIA::PBDTetMeshFEM::solveMaterialConstraint()
{
}

void GAIA::PBDTetMeshFEM::solveBoundaryConstraint()
{
}

void GAIA::PBDTetMeshFEM::solveInversionConstraint()
{
}

GAIA::FloatingType GAIA::PBDTetMeshFEM::evaluateMeritEnergy()
{
	return evaluateInertiaEnergy() + evaluateElasticEnergy();
}

GAIA::FloatingType GAIA::PBDTetMeshFEM::evaluateInertiaEnergy()
{
	float energyInertia = 0.f;

	for (size_t iV = 0; iV < numVertices(); iV++)
	{
		energyInertia += vertexMass(iV) * (mVertPos.col(iV) - inertia.col(iV)).squaredNorm() / (2 * dt * dt);
	}
	return energyInertia;
}



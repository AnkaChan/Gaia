#include "ColliderTetMeshBase.h"

void GAIA::ColliderTetMeshBase::initilizeSurfaceTriMesh()
{
	pSurfaceTriMesh = std::make_shared<TriMeshFEM>();
	pSurfaceTriMesh->resizePositions(numVertices());
	pSurfaceTriMesh->positions() = positions();	
	pSurfaceTriMesh->facePos = pTopology->surfaceFacesTetMeshVIds;

}

void GAIA::ColliderTetMeshBase::updateSurfaceTriMesh()
{
	pSurfaceTriMesh->positions() = positions();
}

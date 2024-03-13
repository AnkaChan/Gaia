#include "VBDPhysicsTest.h"

using namespace GAIA;

__global__ void printGPUCollisionData_kernel(int32_t numActiveCollisions, int32_t *activeCollisionsEachParallelGroup, VBDPhysicsDataGPU* pPhysics) {
    printf("----------------\nnumActiveCollisions: %d\n", numActiveCollisions);
	for (int i = 0; i < numActiveCollisions; i++) {
		int32_t meshId = activeCollisionsEachParallelGroup[i * 2];
		int32_t vertexId = activeCollisionsEachParallelGroup[i * 2 + 1];
		printf("mesh: %d, vertex: %d\n", meshId, vertexId);
		printf("    collision info:\n");
		CollisionDataGPU& collisionData = pPhysics->tetMeshes[meshId]->getCollisionData(vertexId);
		printf("    collide with (%d, %d, %d) from mesh %d:\n", 
			collisionData.closestSurfaceFaceVIds[0], collisionData.closestSurfaceFaceVIds[1], collisionData.closestSurfaceFaceVIds[2], 
			collisionData.intersectedTMeshId);

	}
}

void GAIA::printGPUCollisionData(int32_t numActiveCollisions, int32_t* activeCollisionsEachParallelGroup, VBDPhysicsDataGPU* pPhysics)
{
	printGPUCollisionData_kernel KERNEL_ARGS2(1, 1) (numActiveCollisions, activeCollisionsEachParallelGroup, pPhysics);
	CHECK_CUDA_ERROR(cudaDeviceSynchronize());
}



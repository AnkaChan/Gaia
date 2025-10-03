#pragma once
#include <MeshFrame/Memory/Array.h>


namespace GAIA
{
	struct TriMeshContactRelation
	{
		IdType meshId;
		IdType collisionResultId;	 // to access the collision result list
		IdType collisionId;          // a collition result may have multiple contacts, this is the contact id
		IdType collisionType;        // 0: v-f; 1: edge edge
		IdType collisionVertexOrder; // 0 ~ 4, for v-f contact, 0~3 is from the face side and 4 is from the vertex side, for edge-edge contact, 0~1 is from the first edge and 2~3 is from the second edge
	};
}
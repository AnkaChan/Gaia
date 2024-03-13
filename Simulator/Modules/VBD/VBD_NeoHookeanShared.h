#pragma once

#include "../TetMesh/TetMeshFEMShared.h"

namespace GAIA {
	struct VBD VBDTetMeshNeoHookeanShared : public TetMeshFEMShared
	{
		virtual void initializeGPUTetMesh(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU);
	};
	
	inline void VBDTetMeshNeoHookeanShared::initializeGPUTetMesh(PBDTetMeshMassSpring* pTetMesh, PBDTetMeshFEMGPUMassSpring* pTetMeshGPU)
	{
	}

}
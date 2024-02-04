#include "../TriMesh/TriMesh.h"
#include "TetMeshFEM.h"

namespace GAIA {

	struct TetMeshSurface : protected TriMeshFEM {
		void initialize(TriMeshFEM* pTriMesh);


	};
}
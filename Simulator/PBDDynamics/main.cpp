#include "TetMesh/TetMeshFEM.h"

#include "PBD/PBDPhysics.h"
#include "PBD/PBDTetMeshNeoHookean.h"
#include "Parser/Parser.h"
#include "Parser/InputHandler.h"

#include "CuMatrix/MatrixOps/CuMatrix.h"

int main(int argc, char** argv) {
	REQUIRE_NUM_INPUTS(3);
	GAIA::CommandParser parser;
	GAIA::PBDPhysics physics;
	parser.parse(argc, argv);

	std::string inModelInputFile = argv[1];
	std::string inParameterFile = argv[2];
	std::string outFolder = argv[3];

	InputHandler<GAIA::PBDPhysics> inputHanlder;
	inputHanlder.handleInput(inModelInputFile, inParameterFile, outFolder, parser, physics);

	physics.initializeGPU();

	if (parser.recoveryStateFile != "")
	{
		physics.recoverFromState(parser.recoveryStateFile);
	}

	if (parser.runOnCPU) {
		physics.simulateGPU_debugOnCPU();
	}
	else
	{
		physics.simulateGPU();
	}

}
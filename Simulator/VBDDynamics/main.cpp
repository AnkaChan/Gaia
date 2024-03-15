
#include "Parser/Parser.h"
#include "Parser/InputHandler.h"
#include "TetMesh/TetMeshFEM.h"

#include "VBD/VBDPhysics.h"

int main(int argc, char** argv) {
	REQUIRE_NUM_INPUTS(3);
	GAIA::CommandParser parser;
	GAIA::VBDPhysics physics;
	parser.parse(argc, argv);

	std::string inModelInputFile = argv[1];
	std::string inParameterFile = argv[2];
	std::string outFolder = argv[3];

	InputHandlerEBD<GAIA::VBDPhysics> inputHanlder;
	inputHanlder.handleInput(inModelInputFile, inParameterFile, outFolder, parser, physics);

	physics.initialize();

	if (parser.recoveryStateFile != "")
	{
		physics.recoverFromState(parser.recoveryStateFile);
	}

	physics.simulate();
}
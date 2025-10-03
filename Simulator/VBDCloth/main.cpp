#include <TetMesh/TetMeshFEM.h>

#include <VBDCloth/VBDClothPhysics.h>
#include "Parser/Parser.h"
#include "Parser/InputHandler.h"

#include "CuMatrix/MatrixOps/CuMatrix.h"

#include "Test.h"

int main(int argc, char** argv) {
	// simulateClothMeshStVK();

	test();

	REQUIRE_NUM_INPUTS(3);
	GAIA::CommandParser parser;
	GAIA::VBDClothSimulationFramework physics;
	parser.parse(argc, argv);

	std::string inModelInputFile = argv[1];
	std::string inParameterFile = argv[2];
	std::string outFolder = argv[3];

	InputHandlerVBDCloth<GAIA::VBDClothSimulationFramework> inputHanlder;
	inputHanlder.handleInput(inModelInputFile, inParameterFile, outFolder, parser, physics);

	physics.initialize();

	if (parser.recoveryStateFile != "")
	{
		physics.recoverFromState(parser.recoveryStateFile);
	}

	physics.simulate();

}
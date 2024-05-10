#include <GraphColoring/Graph.h>

#include <GraphColoring/TetMeshTetGraph.h>
#include <GraphColoring/TetMeshEdgeGraph.h>
#include <GraphColoring/TetMeshVertexGraph.h>

#include <GraphColoring/TriMeshVertexGraph.h>

#include <GraphColoring/mcs.h>
#include <GraphColoring/gready.h>

#include <Parser/Parser.h>
#include <MeshFrame/Utility/Parser.h>
#include <MeshFrame/Utility/IO.h>

using namespace GAIA::GraphColoring;

struct GraphColoringParams {
	bool edgeGraph = false;
	bool tGraph = false;
	bool vertexGraph = false;
	bool edgeEnergy = true;
	bool balanceGraphColoring = true;
	float goalMaxMinRatio = 1.05f;

	std::string algorithm = "mcs";

	GraphColoringParams() :
		options("GraphColoring", "Color trimeshes or tetmeshes.")
	{
		options.add_options()
			("e,edgeGraph", "Convert edge constraints to graph, works for trimesh and tetmesh.", cxxopts::value<bool>(edgeGraph))
			("t,tGraph", "Convert tet/triangle based constraints to graph.", cxxopts::value<bool>(tGraph))
			("v,vertexGraph", "Convert vertex connections to graph.", cxxopts::value<bool>(vertexGraph))
			("b,balance", "Whether to balance the size of each color.", cxxopts::value<bool>(balanceGraphColoring))
			("r,ratio", "The goal max/min ratio of the color sizes.", cxxopts::value<float>(goalMaxMinRatio))
			("E,edgeEnergy", "Whether to added edge based bending energy to trimesh vertex graph.", cxxopts::value<bool>(edgeEnergy))
			("a,algorithm", "The name of the coloring algorithm, can be mcs or greedy.", cxxopts::value<std::string>(algorithm))
			;
	}

	void parse(int argc, char** argv) {
		try
		{
			auto result = options.parse(argc, argv);
		}
		catch (const cxxopts::OptionException& e)
		{
			std::cout << "error parsing options: " << e.what() << std::endl;
			std::cout << options.help();
			exit(1);
		}
	}

	cxxopts::Options options;

};

GraphColor::SharedPtr getColoringAlgoritm(const Graph& graph, const GraphColoringParams& params) {

	GraphColor::SharedPtr pColoringAlg = nullptr;

	if (params.algorithm == "mcs")
	{
		pColoringAlg = std::make_shared<Mcs>(graph);
	}
	else if(params.algorithm == "greedy")
	{
		pColoringAlg = std::make_shared<OrderedGreedy>(graph);
	}
	else
	{
		std::cout << "Unsupported algorithm: " << params.algorithm << std::endl;
		exit(1);

	}

	return pColoringAlg;
}

std::shared_ptr<Graph> loadMesh(const std::string& inModelInputFile, const GraphColoringParams& config,
	TMeshStaticF::SharedPtr pTM, TriMeshStaticF::SharedPtr pMesh)
{
	int meshType = -1; // 0: tet mesh; 1: triangular mesh
	MF::IO::FileParts fp(inModelInputFile);

	if (fp.ext == ".t")
	{
		meshType = 0;
		pTM->load_t(inModelInputFile.c_str());
	}
	else if (fp.ext == ".obj")
	{
		meshType = 1;
		pMesh->read_obj(inModelInputFile.c_str());
	}
	else if (fp.ext == ".ply")
	{
		meshType = 1;
		pMesh->read_ply(inModelInputFile.c_str());
	}
	else
	{
		std::cerr << "Error! Unrecognized input extension name: " << fp.ext << std::endl;
		exit(-1);
	}

	std::shared_ptr<Graph> pGraph = nullptr;

	if (config.edgeGraph) {
		std::cout << "Generating edge graph.\n";
		if (meshType == 0)
		{
			pGraph = std::make_shared<TetMeshEdgeGraph>();
			pGraph->fromMesh(pTM.get());
		}
		else if (meshType == 1)
		{
			std::cerr << "Error! Edge graph for triangular mesh is not implemented yet: " << fp.ext << std::endl;
			exit(-1);
		}

	}
	else if (config.tGraph) {
		std::cout << "Generating tet graph.\n";
		if (meshType == 0)
		{
			pGraph = std::make_shared<TetMeshTetGraph>();
			pGraph->fromMesh(pTM.get());
		}
		else if (meshType == 1)
		{
			std::cerr << "Error! Edge graph for triangular mesh is not implemented yet: " << fp.ext << std::endl;
			exit(-1);
		}
	}
	else if (config.vertexGraph)
	{
		if (meshType == 0)
		{
			pGraph = std::make_shared<TetMeshVertexGraph>();
			pGraph->fromMesh(pTM.get());
		}
		else if (meshType) {
			pGraph = std::make_shared<TriMeshVertexGraph>(config.edgeEnergy);
			pGraph->fromMesh(pMesh.get());
		}
	}
	else {
		std::cout << "Please specify which kind of graph you want to generate (edge/tet/vertex).\n";

	}

	return pGraph;
}

int main(int argc, char** argv) {

	REQUIRE_NUM_INPUTS(2);
	GraphColoringParams config;
	config.parse(argc, argv);

	std::string inModelInputFile = argv[1];
	std::string outGraphColoringFile = argv[2];

	TMeshStaticF::SharedPtr pTM = std::make_shared<TMeshStaticF>();
	TriMeshStaticF::SharedPtr pMesh = std::make_shared<TriMeshStaticF>();

	std::shared_ptr<Graph> pGraph = loadMesh(inModelInputFile, config, pTM, pMesh);

	GraphColor::SharedPtr pColoringAlg = getColoringAlgoritm(*pGraph, config);
	
	pColoringAlg->color();

	if (!pColoringAlg->is_valid()) {
		std::cerr << "Graph coloring is invalid" << std::endl;
		return -1;
	}

	pColoringAlg->convertToColoredCategories();

	if (config.balanceGraphColoring)
	{
		pColoringAlg->balanceColoredCategories(config.goalMaxMinRatio);
		if (!pColoringAlg->is_valid()) {
			std::cerr << "Error! Graph coloring is invalid" << std::endl;
			return -1;
		}
	}

	if (outGraphColoringFile.size() != 0)
	{
		pColoringAlg->saveColoringCategories(outGraphColoringFile);
	}

	return 0;
}
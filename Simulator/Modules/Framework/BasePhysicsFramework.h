#pragma once

#include "../Parameters/PhysicsParameters.h"
#include "../Timer/RunningTimeStatistics.h"
#include "../Types/Types.h"
#include "../Materials/Materials.h"
#include "../TetMesh/TetMeshFEM.h"

#include "../Utility/Logger.h"
#include <MeshFrame/Utility/Str.h>
#include <MeshFrame/Utility/IO.h>

namespace GAIA {

	struct CollisionStatistics;
	struct DiscreteCollisionDetector;
	struct ContinuousCollisionDetector;
	struct CollisionDetectionParamters;

	struct ObjectParamsList : public MF::BaseJsonConfig {
		std::vector<ObjectParams::SharedPtr> objectParams;

		size_t size();
		ObjectParams& getObjectParam(int iObj);
		template <typename ObjParamType>
		ObjParamType& getObjectParamAs(int iObj);

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);

		virtual ObjectParams::SharedPtr createObjectParam(const std::string& materialName) = 0;
	};

	struct BasePhysicFramework {
		BasePhysicFramework() {
			simulatorName = "BaseVolumetric";
		}
		typedef std::shared_ptr<BasePhysicFramework> BaseSharedPtr;
		typedef BasePhysicFramework* BasePtr;

		typedef std::shared_ptr<BasePhysicFramework> SharedPtr;
		typedef BasePhysicFramework* Ptr;

		std::vector<std::shared_ptr<TetMeshFEM>> basetetMeshes;
		std::shared_ptr<ObjectParamsList> objectParamsList;
		std::shared_ptr<BasePhysicsParams> basePhysicsParams;
		std::shared_ptr<CollisionDetectionParamters> baseCollisionParams;
		std::shared_ptr<RunningTimeStatistics> baseTimeStatistics;

		std::shared_ptr<DiscreteCollisionDetector> pDCD;
		std::shared_ptr<ContinuousCollisionDetector> pCCD;

		void updateWorldBox();
		virtual void loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder);
		virtual void parseRunningParameters(nlohmann::json& inModelParams, nlohmann::json& inPhysicsParams);

		virtual void timeStepEqualize();
		virtual void initialize();
		virtual void initializeCollisionDetector();
		virtual void disableModelsLatterToAppear();
		virtual void enableModels();

		virtual void recoverFromState(std::string& stateFile);

		virtual void simulate();
		virtual void runStep() = 0;

		virtual std::string getDebugFolder();

		// overridden functions: must be overriden to create certain materials
		virtual TetMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, TetMeshMF::SharedPtr pTMeshMF, BasePhysicsParams::SharedPtr physicsParaemters) = 0;
		virtual std::shared_ptr<ObjectParamsList> createObjectParamsList() = 0;
		virtual std::shared_ptr<BasePhysicsParams> createPhysicsParams() = 0;
		virtual std::shared_ptr<CollisionDetectionParamters> createCollisionParams();
		virtual std::shared_ptr<RunningTimeStatistics> createRunningTimeStatistics() = 0;

		virtual void setUpOutputFolders(std::string outFolder);
		virtual void writeOutputs(std::string outFolder, int frameId);
		virtual void saveExperimentParameters(const std::string& paramsOutOutPath, int indent = 2);

		virtual bool writeSimulationParameters(nlohmann::json& outPhysicsParams);

		virtual void debugOperation(int debugLvl, std::function<void()> ops);
		virtual void debugPrint(int debugLvl, std::string info);
		virtual void saveDebugState(const std::string customName, bool saveMesh=false, const std::string outfolder = "");

		size_t numMeshes() { return basetetMeshes.size(); };


		// time and frames
		int frameId = 0;
		int substep;
		int iIter;
		FloatingType curTime = 0;

		// inputs and outputs
		std::string inputModelListFile;
		std::string inputParamFile;
		std::string outputFolder;

		std::string simulatorName;

		nlohmann::json modelJsonParams;
		nlohmann::json physicsJsonParams;

		size_t numAllVertices;
		size_t numAllTets;
		size_t numAllEdges;

#ifdef KEEP_MESHFRAME_MESHES
		std::vector<std::shared_ptr<TetMeshMF>> tMeshesMF;
#endif // KEEP_MESHFRAME_MESHES
	};

	template<typename ObjParamType>
	inline ObjParamType& ObjectParamsList::getObjectParamAs(int iObj)
	{
		return *(ObjParamType*)(objectParams[iObj].get());
	}

	struct PhysicsStateMesh : MF::BaseJsonConfig
	{
		std::vector<std::array<double, 3>> velocities;
		std::vector<std::array<double, 3>> position;
		bool fromJson(nlohmann::json& j) {
			EXTRACT_FROM_JSON(j, velocities);
			EXTRACT_FROM_JSON(j, position);
			return true;
		}

		bool toJson(nlohmann::json& j) {
			PUT_TO_JSON(j, velocities);
			PUT_TO_JSON(j, position);
			return true;
		}

	};

	struct PhysicsState : MF::BaseJsonConfig
	{
		std::vector<PhysicsStateMesh> meshesState;
		int frameId = -1;
		double curTime = 0;

		void fromPhysics(const BasePhysicFramework& physics, int inFrameId) {
			frameId = inFrameId;
			curTime = physics.curTime;
			meshesState.clear();
			for (size_t iMesh = 0; iMesh < physics.basetetMeshes.size(); iMesh++)
			{
				meshesState.emplace_back();
				TetMeshFEM::SharedPtr pTM = physics.basetetMeshes[iMesh];

                for (size_t iP = 0; iP < pTM->numVertices(); iP++)
                {
                    std::array<double, 3> velArr = {
                        pTM->mVelocity(0, iP),  pTM->mVelocity(1, iP), pTM->mVelocity(2, iP)
                    };
                    meshesState.back().velocities.push_back(velArr);

					std::array<double, 3> ptArr = {
						pTM->mVertPos(0, iP),  pTM->mVertPos(1, iP), pTM->mVertPos(2, iP)
					};
					meshesState.back().position.push_back(ptArr);
				}
			}
		}

		void initializePhysics(BasePhysicFramework& physics, int& inFrameId) {
			inFrameId = frameId;
			physics.curTime = curTime;
			for (size_t iMesh = 0; iMesh < physics.basetetMeshes.size(); iMesh++)
			{
				TetMeshFEM::SharedPtr pTM = physics.basetetMeshes[iMesh];
				if (iMesh >= meshesState.size())
				{
					break;
				}

				for (size_t iP = 0; iP < pTM->numVertices(); iP++)
				{
					pTM->mVelocity(0, iP) = meshesState[iMesh].velocities[iP][0];
					pTM->mVelocity(1, iP) = meshesState[iMesh].velocities[iP][1];
					pTM->mVelocity(2, iP) = meshesState[iMesh].velocities[iP][2];

					pTM->mVertPos(0, iP) = meshesState[iMesh].position[iP][0];
					pTM->mVertPos(1, iP) = meshesState[iMesh].position[iP][1];
					pTM->mVertPos(2, iP) = meshesState[iMesh].position[iP][2];
				}
			}
		}

		bool fromJson(nlohmann::json& j) {
			EXTRACT_FROM_JSON(j, frameId);
			EXTRACT_FROM_JSON(j, curTime);

			for (nlohmann::json& mesh : MF::tryGetJson(j, "meshesState"))
			{
				meshesState.emplace_back();
				meshesState.back().fromJson(mesh);
			}

			return true;
		}

		bool toJson(nlohmann::json& j) {

			nlohmann::json meshesStateJson;
			for (PhysicsStateMesh& mesh : meshesState)
			{
				meshesStateJson.emplace_back();
				mesh.toJson(meshesStateJson.back());
			}

			j["meshesState"] = meshesStateJson;
			j["frameId"] = frameId;
			j["curTime"] = curTime;
			j["format"] = "text"; // ["text", "binary"]

			return true;
		}

		bool writeToJsonFile(std::string filePath, BasePhysicFramework& physics, int& inFrameId, bool binary = true) {
			if (!binary) {
				fromPhysics(physics, inFrameId);
				return BaseJsonConfig::writeToJsonFile(filePath);
			}
			nlohmann::json j;
			j["frameId"] = inFrameId;
			j["curTime"] = physics.curTime;
			j["format"] = "binary"; // ["text", "binary"]
			int nMeshes = physics.basetetMeshes.size();
			j["nMeshes"] = nMeshes;
			nlohmann::json meshesStateJson;
			std::vector<int> nVerts;
			auto filePart = MF::IO::FileParts(filePath);
			std::string basename = filePart.name;
			std::string dirPath = filePart.path;
			std::string velocities_file = basename + ".vel";
			std::string position_file = basename + ".pos";
			meshesStateJson["velocities_file"] = velocities_file;
			meshesStateJson["position_file"] = position_file;
			try {
				// create the file and use binary and append mode
				std::ofstream velFile(dirPath + "/" + velocities_file, std::ios::binary);
				std::ofstream posFile(dirPath + "/" + position_file, std::ios::binary);
				for (int iMesh = 0; iMesh < nMeshes; iMesh++)
				{
					TetMeshFEM::SharedPtr pTM = physics.basetetMeshes[iMesh];
					nVerts.push_back(pTM->numVertices());
					// append to the binary files

					velFile.write((char*)pTM->velocities().data(), sizeof(FloatingType) * pTM->velocities().size());
					posFile.write((char*)pTM->vertices().data(), sizeof(FloatingType) * pTM->vertices().size());
				}
			}
			catch (std::exception& e) {
				std::cout << "Error writing to file: " << e.what() << std::endl;
				return false;
			}
			meshesStateJson["nVerts"] = nVerts;
			j["meshesState"] = meshesStateJson;
			bool retVal = MF::saveJson(filePath, j, 4);
			return retVal;
		}

		bool loadFromJsonFile(std::string filePath, BasePhysicFramework& physics, int& inFrameId) {
			nlohmann::json j;
			MF::loadJson(filePath, j);
			std::string format = "text";
			EXTRACT_FROM_JSON(j, format);

			if (format == "binary") {
				EXTRACT_FROM_JSON(j, frameId);
				EXTRACT_FROM_JSON(j, curTime);
				inFrameId = frameId;
				physics.curTime = curTime;
				int nMeshes{};
				EXTRACT_FROM_JSON(j, nMeshes);
				nMeshes = std::min(nMeshes, (int)physics.basetetMeshes.size());

				nlohmann::json meshesStateJson = MF::tryGetJson(j, "meshesState");
				std::vector<int> nVerts{};
				EXTRACT_FROM_JSON(meshesStateJson, nVerts);
				assert(nVerts.size() >= nMeshes);
				std::string velocities_file;
				std::string position_file;
				EXTRACT_FROM_JSON(meshesStateJson, velocities_file);
				EXTRACT_FROM_JSON(meshesStateJson, position_file);

				auto filePart = MF::IO::FileParts(filePath);
				auto dir_path = filePart.path;
				// load the binary files
				try {
					std::ifstream velFile(dir_path + "/" + velocities_file, std::ios::binary);
					std::ifstream posFile(dir_path + "/" + position_file, std::ios::binary);
					for (int iMesh = 0; iMesh < nMeshes; iMesh++) {
						TetMeshFEM::SharedPtr pTM = physics.basetetMeshes[iMesh];
						int nVert = nVerts[iMesh];
						if (nVert !=pTM->numVertices())
						{
							std::cout << "Error! Number of vertices from recovery state doesn't match! "<< std::endl;
							assert(false);
							getchar();
						}
						velFile.read((char*)pTM->velocities().data(), sizeof(FloatingType) * pTM->velocities().size());
						posFile.read((char*)pTM->vertices().data(), sizeof(FloatingType) * pTM->vertices().size());
					}
				}
				catch (std::exception& e) {
					std::cout << "Error reading from file: " << e.what() << std::endl;
					return false;
				}
				return true;
			}
			else if (format == "text") {
				bool retVal = BaseJsonConfig::loadFromJsonFile(filePath);
				if (retVal) {
					initializePhysics(physics, inFrameId);
				}
				return retVal;

			}

		}
	};
}
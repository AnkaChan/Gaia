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
		std::vector<std::array<FloatingType, 3>> velocities;
		std::vector<std::array<FloatingType, 3>> position;
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

	struct PhysicsState : public MF::BaseJsonConfig
	{
		std::vector<PhysicsStateMesh> meshesState;
		int frameId = -1;
		double curTime = 0;
		bool binary = false;

		void fromPhysics(const BasePhysicFramework& physics) {
			frameId = physics.frameId;
			curTime = physics.curTime;
			meshesState.clear();
			for (size_t iMesh = 0; iMesh < physics.basetetMeshes.size(); iMesh++)
			{
				meshesState.emplace_back();
				TetMeshFEM::SharedPtr pTM = physics.basetetMeshes[iMesh];

                for (size_t iP = 0; iP < pTM->numVertices(); iP++)
                {
                    std::array<FloatingType, 3> velArr = {
                        pTM->mVelocity(0, iP),  pTM->mVelocity(1, iP), pTM->mVelocity(2, iP)
                    };
                    meshesState.back().velocities.push_back(velArr);

					std::array<FloatingType, 3> ptArr = {
						pTM->mVertPos(0, iP),  pTM->mVertPos(1, iP), pTM->mVertPos(2, iP)
					};
					meshesState.back().position.push_back(ptArr);
				}
			}
		}

		void initializePhysics(BasePhysicFramework& physics) {
			physics.frameId = frameId;
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

			std::string format = "text";
			EXTRACT_FROM_JSON(j, format);

			if (format == "binary") {
				assert(userData != nullptr);
				const std::string & filePath = *(std::string *)(userData);

				int nMeshes{};
				EXTRACT_FROM_JSON(j, nMeshes);

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
						int nVert = nVerts[iMesh];
						meshesState.emplace_back();

						PhysicsStateMesh & meshState = meshesState.back();
						meshState.velocities.resize(nVert);
						meshState.position.resize(nVert);

						velFile.read((char*)meshState.velocities.data(), sizeof(FloatingType) * nVert * 3);
						posFile.read((char*)meshState.position.data(), sizeof(FloatingType) * nVert * 3);
					}
				}
				catch (std::exception& e) {
					std::cout << "Error reading from file: " << e.what() << std::endl;
					return false;
				}
				return true;
			}
			else if (format == "text") {
				for (nlohmann::json& mesh : MF::tryGetJson(j, "meshesState"))
				{
					meshesState.emplace_back();
					meshesState.back().fromJson(mesh);
				}
				
				return true;

			}

			return true;
		}

		bool toJson(nlohmann::json& j) {
			j["frameId"] = frameId;
			j["curTime"] = curTime;

			if (binary)
			{
				assert(userData != nullptr);
				const std::string& filePath = *(std::string*)(userData);

				j["format"] = "binary"; // ["text", "binary"]
				size_t nMeshes = meshesState.size();
				j["nMeshes"] = meshesState.size();
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
						int nVert = meshesState[iMesh].velocities.size();
						nVerts.push_back(nVert);
						const auto& meshState = meshesState[iMesh];
						velFile.write((char*)meshState.velocities.data(), sizeof(FloatingType) * nVert * 3);
						posFile.write((char*)meshState.position.data(), sizeof(FloatingType) * nVert * 3);
					}
				}
				catch (std::exception& e) {
					std::cout << "Error writing to file: " << e.what() << std::endl;
					return false;
				}
				meshesStateJson["nVerts"] = nVerts;
				j["meshesState"] = meshesStateJson;
				return true;
			}
			else
			{
				j["format"] = "text"; // ["text", "binary"]
				nlohmann::json meshesStateJson;
				for (PhysicsStateMesh& mesh : meshesState)
				{
					meshesStateJson.emplace_back();
					mesh.toJson(meshesStateJson.back());
					j["meshesState"] = meshesStateJson;
				}
			}


			return true;
		}

	};
}
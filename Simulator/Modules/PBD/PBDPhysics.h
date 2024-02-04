#pragma once

#include "../Types/Types.h"
#include "PBDTetMeshFEM.h"
#include "PBDTetMeshFEMGPU.h"

#include "PBDPhysicsParameters.h"
#include "../CollisionDetector/CollisionDetertionParameters.h"
#include "Timer/Timer.h"
#include "Timer/RunningTimeStatistics.h"
#include "Deformer.h"

#include "PBDTetMeshNeoHookean.h"
#include "PBDTetMeshMassSpring.h"

#include "cuda_runtime.h"

// #define DO_COLLISION_STATISTICS 

namespace GAIA {
    struct DiscreteCollisionDetector;
    struct ContinuousCollisionDetector;
    struct VolumetricCollisionDetector;

	struct ObjectsParamList : public MF::BaseJsonConfig {
		std::vector<ObjectParamsPBD::SharedPtr> objectParams;

		size_t size();
		ObjectParamsPBD& getObjectParam(int iObj);

		virtual bool fromJson(nlohmann::json& objectParam);
		virtual bool toJson(nlohmann::json& objectParam);
	};

	struct PBDPhysicsAllParameters : public MF::BaseJsonConfig {
		PBDPhysicParameters physicsParams;
		CollisionDetectionParamters collisionParams;

		virtual bool fromJson(nlohmann::json& physicsJsonParams);
		virtual bool toJson(nlohmann::json& physicsJsonParams);
	};


	/*
	* = PBDPhysics + SoftBodyManager + GeoManager
	*/
	struct PBDPhysics
    {
        typedef std::shared_ptr<PBDTetMeshFEM> BaseMaterialPtr;
        typedef std::shared_ptr<PBDTetMeshFEMShared_NeoHookean> NeoHookeanMaterialPtr;
        typedef std::shared_ptr<PBDTetMeshFEMShared_MassSpring> MassSpringMaterialPtr;

		// inModelInputFile: a json file that contains input model files and model-specific configurations
		// inParameterFile:  a json file contains global parameters that are used by collision handler and ccd
		void loadRunningparameters(std::string inModelInputFile, std::string inParameterFile, std::string outFolder);

        void saveExperimentParameters(const std::string& paramsOutOutPath);

		bool initializeGPU();

        void updateGPUMeshes();

		void syncAllToGPU(bool sync = true);
        void syncAllToCPU(bool sync = true);
        void syncAllToGPUVertPosOnly(bool sync = true);
        void syncAllToCPUVertPosOnly(bool sync = true);
        void syncAllToInvertedSignOnly(bool sync = true);

        void recoverFromState(std::string& stateFile);

		void simulateGPU();
		void simulateGPU_debugOnCPU();

		void runStepGPU();
		void runStepGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU);
		
		void materialBoundarySolveGPU();
		void materialBoundarySolveGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU);

        void InversionSolveGPU();

		void BoundarySolve();
		void BoundarySolveSolveGPU_debugOnCPU(std::vector<PBDTetMeshFEMGPU*>& gpuTMeshesOnCPU);

        void collisionSolve();
        void depthBasedCollisionSolve();
        void volumeBasedCollisionSolve();

        void initialCollisionsDetection();
        void updateDCD_CPU(bool rebuild);
        void updateCCD_CPU(bool rebuild);

        void changeNumSubSteps(int newNumSubSteps);

        void solveCollisionConstraintsCPU();
        void solveCollisionConstraintsOneVertCPU(int32_t curVertID, int32_t curMeshID, 
            CollisionDetectionResult& colResult, int32_t iIntersection);
        void applyCollisionConstraintsCPU(const Vec3I& closestFaceVIds, int32_t curvertId, int32_t curMeshID, int32_t collidedMeshID,
            FloatingType C, FloatingType penetrationDepth, FloatingType dt, Vec3& barycentric, const Vec12& gradients);

		void updateVelocities();

        void updateWorldBox();

		void writeOutputs(std::string outFolder, int frameId);

		inline ObjectParamsPBD& getObjectParam(int iObj) {
			return objectParamsList.getObjectParam(iObj);
		}

		inline GAIA::PBDPhysicParameters& GAIA::PBDPhysics::physicsParams()
		{
			return physicsAllParams.physicsParams;
		}

        inline GAIA::CollisionDetectionParamters& GAIA::PBDPhysics::collisionParams()
        {
            return physicsAllParams.collisionParams;
        }

		bool initializeFromState(std::string& stateFile);

		void setUpOutputFolders(std::string outFolder);

		void applyInitialGuesses();

		void evaluateObjectActiveness();

        void applyDeformers();
        void applyConstraints();
        void applyPostSolvingDeformers();

        void solveInvertedTets();
        void solveInversionConstraint(PBDTetMeshFEM * pTM, int iTet);

        // evaluation and statistics
        void evaluateConvergence(int iStep, int iIteration);

        std::vector<DeformerPtr>& getDeformers();
        std::vector<DeformerPtr>& getPostSolvingDeformers();


		// tetmesh and corresponding materials
		std::vector<std::shared_ptr<PBDTetMeshFEM>> tMeshes;
#ifdef KEEP_MESHFRAME_MESHES
        std::vector<std::shared_ptr<TetMeshMF>> tMeshesMF;
#endif // KEEP_MESHFRAME_MESHES

        std::shared_ptr<DiscreteCollisionDetector> pDCD;
        std::shared_ptr<ContinuousCollisionDetector> pCCD;
        std::shared_ptr<VolumetricCollisionDetector> pVolCD;


        std::vector<PBDTetMeshFEMGPU*> GPUTMeshes;

		PBDPhysicsAllParameters physicsAllParams;

		ObjectsParamList objectParamsList;
		// time and frames
		int frameId = 0;
        int substep;
        int iIter;
        int collisionDetectionCounter = 0;

		FloatingType curTime = 0;
		FloatingType dt;
		RunningTimeStatistics timeStatistics;

		// std::vector<DeformerPtr> deformers;
		// std::vector<DeformerPtr> postSolvingDeformers;

		std::string inputModelListFile;
		std::string inputParamFile;
		std::string outputFolder;
        //
        std::vector<std::vector<IdType>> initialPenetratedVertIds;
        std::vector<std::vector<IdType>> initialPenetrationFreeVertIds;

		std::vector<cudaStream_t> cudaStreams;

        std::vector<CollisionDetectionResult> positiveCollisionDetectionResults;

        std::vector<std::vector<CollisionDetectionResult>> ccdResultsAll;
        std::vector<std::vector<CollisionDetectionResult>> dcdResultsAll;
        std::vector<DeformerPtr> deformers;
        std::vector<DeformerPtr> postSolvingDeformers;

#ifdef DO_COLLISION_STATISTICS 
        CollisionStatistics collisionStatistics;

#endif //  DO_COLLISION_STATISTICS 
    


	};

    struct PBDPhysicsStateMesh : MF::BaseJsonConfig
    {
        std::vector<std::array<double, 3>> velocity;
        std::vector<std::array<double, 3>> position;
        bool fromJson(nlohmann::json& j) {
            EXTRACT_FROM_JSON(j, velocity);
            EXTRACT_FROM_JSON(j, position);
            return true;
        }

        bool toJson(nlohmann::json& j) {
            PUT_TO_JSON(j, velocity);
            PUT_TO_JSON(j, position);
            return true;
        }

    };

    struct PBDPhysicsState : MF::BaseJsonConfig
    {
        std::vector<PBDPhysicsStateMesh> meshesState;
        int frameId = -1;
        double curTime = 0;

        void fromPhysics(const PBDPhysics& physics, int inFrameId) {
            frameId = inFrameId;
            curTime = physics.curTime;
            meshesState.clear();
            for (size_t iMesh = 0; iMesh < physics.tMeshes.size(); iMesh++)
            {
                meshesState.emplace_back();
                PBDTetMeshFEM::SharedPtr pTM = physics.tMeshes[iMesh];

                for (size_t iP = 0; iP < pTM->mVelocity.cols(); iP++)
                {
                    std::array<double, 3> velArr = {
                        pTM->mVelocity(0, iP),  pTM->mVelocity(1, iP), pTM->mVelocity(2, iP)
                    };
                    meshesState.back().velocity.push_back(velArr);

                    std::array<double, 3> ptArr = {
                        pTM->mVertPos(0, iP),  pTM->mVertPos(1, iP), pTM->mVertPos(2, iP)
                    };
                    meshesState.back().position.push_back(ptArr);
                }
            }
        }

        void initializeSoftBodyManager(PBDPhysics& physics, int& inFrameId) {
            inFrameId = frameId;
            physics.curTime = curTime;
            for (size_t iMesh = 0; iMesh < physics.tMeshes.size(); iMesh++)
            {
                PBDTetMeshFEM::SharedPtr pTM = physics.tMeshes[iMesh];
                if (iMesh >= meshesState.size())
                {
                    break;
                }

                for (size_t iP = 0; iP < pTM->numVertices(); iP++)
                {
                    pTM->mVelocity(0, iP) = meshesState[iMesh].velocity[iP][0];
                    pTM->mVelocity(1, iP) = meshesState[iMesh].velocity[iP][1];
                    pTM->mVelocity(2, iP) = meshesState[iMesh].velocity[iP][2];

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
            for (PBDPhysicsStateMesh& mesh : meshesState)
            {
                meshesStateJson.emplace_back();
                mesh.toJson(meshesStateJson.back());
            }

            j["meshesState"] = meshesStateJson;
            j["frameId"] = frameId;
            j["curTime"] = curTime;

            return true;
        }
    };


	inline GAIA::ObjectParamsPBD& GAIA::ObjectsParamList::getObjectParam(int iObj)
	{
		return *objectParams[iObj];
	}
}
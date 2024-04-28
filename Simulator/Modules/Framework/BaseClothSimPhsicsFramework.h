#pragma once

#include "BasePhysicsFramework.h"
#include "../TriMesh/TriMesh.h"
#include "../SpatialQuery/DynamicCollider.h"
#include "../SpatialQuery/ColiiderTriMeshBase.h"
#include <CollisionDetector/ClothContactDetector.h>

namespace GAIA {
	struct BaseClothPhsicsFramework : public BasePhysicFramework
	{
		typedef std::shared_ptr<BaseClothPhsicsFramework> SharedPtr;
		typedef BaseClothPhsicsFramework* Ptr;

		virtual TriMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, BasePhysicsParams::SharedPtr physicsParaemters, BaseClothPhsicsFramework* pPhysics) = 0;
		// deprecate the volumetric simulation framework
		virtual TetMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, TetMeshMF::SharedPtr pTMeshMF, BasePhysicsParams::SharedPtr physicsParaemters);

        virtual void parseRunningParameters() override;

		virtual void initialize();
        // collider includes meshes not from simulation, such as  
        virtual void initializeCollider();
        // return a pointer to the base class of the collider mesh
        ColliderTrimeshBase::SharedPtr createColliderMesh(nlohmann::json& colliderMeshJsonParams);
        virtual void initializeViewer();

        virtual void updateCollider();

		virtual void writeOutputs(std::string outFolder, int frameId);
		virtual void runStep()=0;
        virtual void recoverFromState(std::string& stateFile);

		virtual size_t numTriMeshes() { return baseTriMeshesForSimulation.size(); };

        size_t numSimulationMeshes() { return baseTriMeshesForSimulation.size(); };

        virtual void setSimulatedMeshToUpToDateStatus(bool updated);

		bool isSimulationMesh(IdType meshId) {
			return meshId < baseTriMeshesForSimulation.size();
		}

        std::vector<TriMeshFEM::SharedPtr> baseTriMeshesForSimulation;
        std::vector<ColliderTrimeshBase::SharedPtr> colliderMeshes;

        // baseTriMeshesForSimulation + colliderMeshes, in that order
        std::vector<TriMeshFEM::SharedPtr> triMeshesAll;
        
        ClothContactDetectorParameters::SharedPtr pClothContactDetectorParameters;
        ClothContactDetector::SharedPtr pClothContactDetector;

        DynamicColliderParameters::SharedPtr pDynamicColliderParameter;
        DynamicCollider::SharedPtr pDynamicCollider;

        nlohmann::json colliderJsonParams;
        nlohmann::json contactDetectorParams;



	};
	
	struct ClothPhysicsState : public PhysicsState
	{
        void fromPhysics(const BasePhysicFramework& physics, int inFrameId) = delete;
		void toPhysics(BasePhysicFramework& physics, int inFrameId) = delete;
        
        void fromPhysics(const BaseClothPhsicsFramework& physics) {
            frameId = physics.frameId;
            curTime = physics.curTime;
            meshesState.clear();
            for (size_t iMesh = 0; iMesh < physics.baseTriMeshesForSimulation.size(); iMesh++)
            {
                meshesState.emplace_back();
                TriMeshFEM::SharedPtr pMesh = physics.baseTriMeshesForSimulation[iMesh];

                for (size_t iP = 0; iP < pMesh->velocities.cols(); iP++)
                {
                    std::array<FloatingType, 3> velArr = {
                        pMesh->velocities(0, iP),  pMesh->velocities(1, iP), pMesh->velocities(2, iP)
                    };
                    meshesState.back().velocities.push_back(velArr);

                    std::array<FloatingType, 3> ptArr = {
                        pMesh->positions()(0, iP),  pMesh->positions()(1, iP), pMesh->positions()(2, iP)
                    };
                    meshesState.back().position.push_back(ptArr);
                }
            }
        }

        void initializePhysics(BaseClothPhsicsFramework& physics) {
            physics.frameId = frameId;
            physics.curTime = curTime;
            for (size_t iMesh = 0; iMesh < physics.baseTriMeshesForSimulation.size(); iMesh++)
            {
                TriMeshFEM::SharedPtr pMesh = physics.baseTriMeshesForSimulation[iMesh];;
                if (iMesh >= meshesState.size())
                {
                    break;
                }

                for (size_t iP = 0; iP < pMesh->numVertices(); iP++)
                {
                    pMesh->velocities(0, iP) = meshesState[iMesh].velocities[iP][0];
                    pMesh->velocities(1, iP) = meshesState[iMesh].velocities[iP][1];
                    pMesh->velocities(2, iP) = meshesState[iMesh].velocities[iP][2];

                    pMesh->positions()(0, iP) = meshesState[iMesh].position[iP][0];
                    pMesh->positions()(1, iP) = meshesState[iMesh].position[iP][1];
                    pMesh->positions()(2, iP) = meshesState[iMesh].position[iP][2];
                }
            }
        }
	};
}

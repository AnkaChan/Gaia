#pragma once

#include "BasePhysicsFramework.h"
#include "../TriMesh/TriMesh.h"

namespace GAIA {
	struct BaseClothPhsicsFramework : public BasePhysicFramework
	{
		typedef std::shared_ptr<BaseClothPhsicsFramework> SharedPtr;
		typedef BaseClothPhsicsFramework* Ptr;

		virtual TriMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, BasePhysicsParams::SharedPtr physicsParaemters, BaseClothPhsicsFramework* pPhysics) = 0;
		// deprecate the volumetric simulation framework
		virtual TetMeshFEM::SharedPtr initializeMaterial(ObjectParams::SharedPtr objParam, TetMeshMF::SharedPtr pTMeshMF, BasePhysicsParams::SharedPtr physicsParaemters) ;

		virtual void initialize();
		virtual void writeOutputs(std::string outFolder, int frameId);
		virtual void runStep()=0;
        virtual void recoverFromState(std::string& stateFile);


		std::vector<std::shared_ptr<TriMeshFEM>> baseTriMeshes;


	};
	
	struct ClothPhysicsState : public PhysicsState
	{
        void fromPhysics(const BasePhysicFramework& physics, int inFrameId) = delete;
		void toPhysics(BasePhysicFramework& physics, int inFrameId) = delete;
        
        void fromPhysics(const BaseClothPhsicsFramework& physics) {
            frameId = physics.frameId;
            curTime = physics.curTime;
            meshesState.clear();
            for (size_t iMesh = 0; iMesh < physics.baseTriMeshes.size(); iMesh++)
            {
                meshesState.emplace_back();
                TriMeshFEM::SharedPtr pMesh = physics.baseTriMeshes[iMesh];

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
            for (size_t iMesh = 0; iMesh < physics.baseTriMeshes.size(); iMesh++)
            {
                TriMeshFEM::SharedPtr pMesh = physics.baseTriMeshes[iMesh];;
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

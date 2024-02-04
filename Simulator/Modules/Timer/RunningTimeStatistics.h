#pragma once
#include <MeshFrame/Utility/Parser.h>" << std::endl
#include <sstream>
namespace GAIA {
    struct ConvergenceStats : MF::BaseJsonConfig {
        ConvergenceStats(int numMeshes)
        {
            energy.resize(numMeshes, 0.f);
            energy_inertia.resize(numMeshes, 0.f);
            energy_elastic.resize(numMeshes, 0.f);
            energy_elastic_StVK.resize(numMeshes, 0.f);
            energy_elastic_bending.resize(numMeshes, 0.f);
            energy_collision.resize(numMeshes, 0.f);
            energy_registration.resize(numMeshes, 0.f);
            avgGradNorm.resize(numMeshes, 0.f);
            energy_collision_body.resize(numMeshes, 0.f);
        }

        // by meshes
        FloatingType energy_allClothes;
        FloatingType avgGradNorm_allClothes;

        std::vector<FloatingType> energy;
        std::vector<FloatingType> energy_inertia;
        std::vector<FloatingType> energy_elastic;
        std::vector<FloatingType> energy_elastic_StVK;
        std::vector<FloatingType> energy_elastic_bending;
        std::vector<FloatingType> energy_collision;
        std::vector<FloatingType> energy_collision_body;
        std::vector<FloatingType> energy_registration;
        std::vector<FloatingType> avgGradNorm;
        
        FloatingType stepSize = -1;
        FloatingType initialStepSize = -1;

        bool fromJson(nlohmann::json& objectParam) {
            return true;
        }

        bool toJson(nlohmann::json& objectParam) {
            PUT_TO_JSON(objectParam, energy);
            PUT_TO_JSON(objectParam, energy_inertia);
            PUT_TO_JSON(objectParam, energy_elastic);
            PUT_TO_JSON(objectParam, energy_elastic_StVK);
            PUT_TO_JSON(objectParam, energy_elastic_bending);
            PUT_TO_JSON(objectParam, energy_collision);
            PUT_TO_JSON(objectParam, energy_registration);
            PUT_TO_JSON(objectParam, avgGradNorm);

            return true;
        }
    };

    struct RunningTimeStatistics : MF::BaseJsonConfig{
        void setToZero() {
            timeCsmpFrame = 0;

            timeCsmpAllSubSteps = 0;

            timeCsmpInitialStep = 0;
            timeCsmpMaterialSolve = 0;
            timeCsmpCollisionSolve = 0;
            timeCsmpInversionSolve = 0;

            timeCsmpInitialCollision = 0;

            timeCsmpUpdatingCollisionInfoDCD = 0;
            timeCsmpUpdatingBVHDCD = 0;
            timeCsmpColDetectDCD = 0;
            timeCsmpShortestPathSearchDCD = 0;

            timeCsmpUpdatingCollisionInfoCCD = 0;
            timeCsmpUpdatingBVHCCD = 0;
            timeCsmpColDetectCCD = 0;
            timeCsmpUpdateVelocity = 0;
            timeCsmpSaveOutputs = 0;

            meritEnergy.clear();
        }

        virtual std::string customString() {
            return std::string();
        }

        virtual std::string getString() {
            std::stringstream ss;
            ss << "---FRAME INFORMATION" << "\n";
            ss << "-----Step Total: " << timeCsmpAllSubSteps << "\n";
            ss << "-----Initial Step: " << timeCsmpInitialStep << "\n";
            ss << "-----Material Solve: " << timeCsmpMaterialSolve << "\n";
            ss << "-----Initial Collision Detection: " << timeCsmpInitialCollision << "\n";
            ss << "-----Inversion Solve: " << timeCsmpInversionSolve << "\n";
            ss << "-----DCD Collision Information Uptate: " << timeCsmpUpdatingCollisionInfoDCD << "\n";
            ss << "---------DCD Uptating BVH: " << timeCsmpUpdatingBVHDCD << "\n";
            ss << "---------DCD Detecting Collision: " << timeCsmpColDetectDCD << "\n";
            ss << "---------Shortest Path Search: " << timeCsmpShortestPathSearchDCD << "\n";
            ss << "-----CCD Collision Information Uptate: " << timeCsmpUpdatingCollisionInfoCCD << "\n";
            ss << "---------CCD Uptating BVH: " << timeCsmpUpdatingBVHCCD << "\n";
            ss << "---------CCD Detecting Collision: " << timeCsmpColDetectCCD << "\n";
            ss << "-----Collision Solve: " << timeCsmpCollisionSolve << "\n";
            ss << "-----Updating Velocity: " << timeCsmpUpdateVelocity << "\n";
            ss << customString();
            ss << "-----Save Outputs: " << timeCsmpSaveOutputs << "\n";

            return ss.str();
        }

        void print() {
            std::cout << getString(); 

        }

        virtual bool toJson(nlohmann::json& j)
        {
            PUT_TO_JSON(j, timeCsmpAllSubSteps);

            PUT_TO_JSON(j, timeCsmpMaterialSolve);

            PUT_TO_JSON(j, timeCsmpInitialCollision);

            PUT_TO_JSON(j, timeCsmpInversionSolve);

            PUT_TO_JSON(j, timeCsmpUpdatingCollisionInfoDCD);
            PUT_TO_JSON(j, timeCsmpUpdatingBVHDCD);
            PUT_TO_JSON(j, timeCsmpColDetectDCD);
            PUT_TO_JSON(j, timeCsmpShortestPathSearchDCD);

            PUT_TO_JSON(j, timeCsmpUpdatingCollisionInfoCCD);
            PUT_TO_JSON(j, timeCsmpUpdatingBVHCCD);
            PUT_TO_JSON(j, timeCsmpColDetectCCD);

            PUT_TO_JSON(j, timeCsmpCollisionSolve);

            PUT_TO_JSON(j, timeCsmpSaveOutputs);

            PUT_TO_JSON(j, meritEnergy);

            return true;
        }

        virtual bool fromJson(nlohmann::json& j)
        {
            return true;
        }
        FloatingType timeCsmpFrame = 0;
        FloatingType timeCsmpAllSubSteps = 0;
        
        FloatingType timeCsmpInitialStep = 0;

        FloatingType timeCsmpMaterialSolve = 0;
        FloatingType timeCsmpCollisionSolve = 0;
        FloatingType timeCsmpInversionSolve = 0;

        FloatingType timeCsmpInitialCollision = 0;

        FloatingType timeCsmpUpdatingCollisionInfoDCD = 0;
        FloatingType timeCsmpUpdatingBVHDCD = 0;
        FloatingType timeCsmpColDetectDCD = 0;
        FloatingType timeCsmpShortestPathSearchDCD = 0;

        FloatingType timeCsmpUpdatingCollisionInfoCVolume = 0;
        FloatingType timeCsmpUpdatingBVHVolume = 0;
        FloatingType timeCsmpIntersectionVolumeTest = 0;

        FloatingType timeCsmpUpdatingCollisionInfoCCD = 0;
        FloatingType timeCsmpUpdatingBVHCCD = 0;
        FloatingType timeCsmpColDetectCCD = 0;

        FloatingType timeCsmpUpdateVelocity = 0;


        FloatingType timeCsmpSaveOutputs;


        // convergence statistics
        // step x iteration x object
        std::vector<std::vector<std::vector<FloatingType>>> meritEnergy;
    };


}
from M01_Parameters import *
from M02_GenRunningParameters import *
from M03_ParametersCreator import *

if __name__ == '__main__':
    genFileName = Path(__file__).stem

    machineName = "MyEnv"

    run = True

    models = []

    for i in range(6):
        model2Shifted = getModelInfoOctopusInflatted(modelExample, materialType="MassSpring")
        model2Shifted["translation"] = [0, 20 - 2.5 * i, 0]
        models.append(model2Shifted)


    parameters = getPhysicsParametersoOctopus(parametersExample)

    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.2
    parameters["PhysicsParams"]["boundaryFrictionStatic"] = 0.22
    parameters["PhysicsParams"]["solveInvertedTets"] = True

    parameters["CollisionParams"]["checkFeasibleRegion"] = True
    parameters["CollisionParams"]["feasibleRegionEpsilon"] = 0.02
    parameters["CollisionParams"]["maxSearchDistanceMultiplier"] = 1.5

    parameters["CollisionParams"]["allowCCD"] = True
    parameters["CollisionParams"]["allowDCD"] = True

    parameters["PhysicsParams"]["collisionDetectionSubSteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 2
    parameters["PhysicsParams"]["numSubsteps"] = 20
    parameters["PhysicsParams"]["numFrames"] = 300

    parameters["PhysicsParams"]["worldBounds"] = [
        [-2.0, 0.0, -2.0],
        [ 2.0, 200.0, 2.0]
    ]
    modelsInfo = {
        "Models":models
    }

    experimentName = "Test_Substep20_DCD_CCD"
    parameters["CollisionParams"]["allowCCD"] = True
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["PhysicsParams"]["iterations"] = 3
    parameters["PhysicsParams"]["numSubsteps"] = 20
    genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters,  runCommand=run)




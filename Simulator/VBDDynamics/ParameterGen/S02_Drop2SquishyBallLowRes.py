from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np

def selectFixedPoints(initialStateFile, heightThres=4.2):
    data = json.load(open(initialStateFile))
    pts = np.array(data["meshesState"][0]["position"])

    selectedIds = np.where(pts[:, 1]>heightThres)

    return selectedIds[0].tolist()

if __name__ == '__main__':
    genFileName = Path(__file__).stem

    selectedIds = []
    machineName = "AnkaPC00"
    binaryFile = machines[machineName]["binaryFile"]

    # run = False
    run = True

    models = []
    numModels = 2

    modelsLow = 8
    modelsHigh = 20
    vDisplacement = 1
    fixedPoints = selectedIds

    for iModel in range(numModels):
        model = getModelInfoSquishyBallLowRes(modelExample)
        model['density'] = 20
        model["fixedPoints"] = fixedPoints

        model["translation"] = [vDisplacement * iModel, modelsLow + iModel * (modelsHigh - modelsLow) / (numModels-1), 0]
        model['scale'] = [0.3, 0.3, 0.3]
        model['initialVelocity'] = [0, -5, 0]
        model["exponentialVelDamping"] = 0.95
        model["constantVelDamping"] = 0.02
        model["frictionDynamic"] = 0.1
        model['dampingHydrostatic'] = 1e-7
        model["dampingDeviatoric"] = 1e-7

        # model["epsV"] = 1
        model["frictionEpsV"] = 0.01

        models.append(model)

    modelsInfo = {
        "Models": models
    }

    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 3000
    parameters["PhysicsParams"]["worldBounds"] = [[-10, 0, -10], [10, 20, 10]]
    recoveryState = None

    experimentName = "Test_2Ball_step10_iter10_0.0ColThickness_GPU"
    parameters["PhysicsParams"]["gravity"] = [0.0, -10.0, 0.0]
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True
    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.1

    parameters["PhysicsParams"]["useGPU"] = True
    # parameters["PhysicsParams"]["useGPU"] = False

    parameters["CollisionParams"]["restPoseCloestPoint"] = True
    parameters["PhysicsParams"]["collisionStiffness"] = 2e7

    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True

    # turn on to output files
    parameters["PhysicsParams"]["saveOutputs"] = False
    parameters["ViewerParams"] ={
        "enableViewer" : True
    }

    experimentName = "Test_2Ball_step10_iter10_0.0ColThickness_GPU_newImpl_run02"

    parameters["PhysicsParams"]["numSubsteps"] = 4
    parameters["PhysicsParams"]["iterations"] = 24

    parameters["PhysicsParams"]["useAccelerator"]=True
    # parameters["PhysicsParams"]["useAccelerator"]=False
    parameters["PhysicsParams"]["acceleratorRho"] = 0.94
    cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=binaryFile,
                                runCommand=run, recoverState=recoveryState)


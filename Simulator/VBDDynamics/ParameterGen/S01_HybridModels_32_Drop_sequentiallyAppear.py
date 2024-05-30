import math

from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np
from itertools import product

def getModel(i):
    model = None

    density = 50
    maxVelocityMagnitude = 25
    exponentialVelDamping = 0.9
    constantVelDamping = 0.1
    initialVelocity = -5
    modelFriction = 0.5
    modelFrictionEpsV = 1

    if i == 0:
        model = getModelInfoTest_Octupus(modelExample)
        model['scale'] = [0.8, 0.8, 0.8]

    elif i == 1:
        model = getModelInfoTest_BunnyLowRes((modelExample))
        model['scale'] = [0.2, 0.2, 0.2]
        model["translation"] = [0, -1, 0]

    elif i == 2:
        model = getModelInfoTeapot(modelExample)
        model['scale'] = [1.5, 1.5, 1.5]
        model['rotation'] = [-3.1415926/2, 0, 0]
        model['miu'] = 1e5
        model['lmbd'] = 1e6
    elif i == 3:
        model = getModelInfoHighArmadiloLowRes(modelExample)
        model['scale'] = [3, 3, 3]

    else:
        assert False

    model['initialVelocity'] = [0, initialVelocity, 0]
    model['density'] = density
    model["exponentialVelDamping"] = exponentialVelDamping
    model["constantVelDamping"] = constantVelDamping
    model["maxVelocityMagnitude"] = maxVelocityMagnitude
    model["frictionDynamic"] = modelFriction
    model["frictionEpsV"] = modelFrictionEpsV

    model["dampingShear"] = 1e-7
    model["dampingVolume"] = 1e-7
    return model

if __name__ == '__main__':
    genFileName = Path(__file__).stem

    machineName = "AnkaPC00"
    binaryFile = machines[machineName]["binaryFile"]

    # run = False
    run = True

    numModels = 32
    # numModels = 4

    initialElevation = 10
    dropFrameGaps = 6
    startFrame = 120
    elevation = 1.

    r = 2.25

    # fixedPoints = selectedIds.tolist()
    fixedPoints = []
    models = []
    np.random.seed(123)

    for iModel in range(numModels):
        randI = np.random.randint(0,4)
        model = getModel(randI)

        model["translation"] = [r * math.cos(iModel * np.pi / 3), model["translation"][1] + initialElevation, r * math.sin(iModel * np.pi / 3)]

        model["frameToAppear"] = startFrame + iModel * dropFrameGaps
        model["initialVelocity"] = [0, -10, 0]
        model["frictionDynamic"] = 0.2
        model["frictionEpsV"] = 1e-2

        models.append(model)

    modelsInfo = {
        "Models":models
    }
    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 900
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["worldBounds"] = [
      [  0,  0.0,    0 ],
      [ 0,   20.0,   0  ]
    ]

    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True

    parameters["ViewerParams"] ={
        "enableViewer" : True
    }
    # turn on to output files
    parameters["PhysicsParams"]["saveOutputs"] = False

    recoveryState = None

    experimentName = "Test_64HybridModels_2steps_40Iters_stiff_Friction"
    parameters["PhysicsParams"]["useGPU"] = True
    parameters["PhysicsParams"]["useAccelerator"]=True
    # parameters["PhysicsParams"]["useAccelerator"]=False
    parameters["PhysicsParams"]["acceleratorRho"] = 0.98

    # parameters["PhysicsParams"]["useGPU"] = False
    parameters["PhysicsParams"]["collisionStiffness"] = 1e5
    parameters["PhysicsParams"]["numSubsteps"] = 2
    parameters["PhysicsParams"]["iterations"] = 30
    parameters["PhysicsParams"]["collisionSolutionType"] = 1
    # recoveryState = r"E:\Data2\VBDSimulation\S09_Experiment_Hybrid64_sequentiallyAppear\Test_64HybridModels_4steps_15Iters_stiff_Friction\RecoveryStates\A00000050.json"
    cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=binaryFile, runCommand=run, recoverState=recoveryState)
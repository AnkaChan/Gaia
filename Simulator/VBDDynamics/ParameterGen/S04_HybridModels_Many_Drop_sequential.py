import math

from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np
from itertools import product

def getModel(i):
    model = None

    density = 50
    maxVelocityMagnitude = 25
    exponentialVelDamping = 0.95
    constantVelDamping = 0.01
    initialVelocity = -5
    frictionDynamic = 0.15

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

    model["dampingHydrostatic"] = 1e-7
    model["dampingDeviatoric"] = 1e-7
    model['frictionDynamic'] = frictionDynamic
    model['initialVelocity'] = [0, initialVelocity, 0]
    model['density'] = density
    model["exponentialVelDamping"] = exponentialVelDamping
    model["constantVelDamping"] = constantVelDamping
    model["maxVelocityMagnitude"] = maxVelocityMagnitude

    return model

def addModelsCuboid(newModels, hStart, w, l, h, vDis=5.6, horizontalDis=1.5, layerStartToRelaseSequentially=10, releaseFrameGap=15):

    xTranslations = np.linspace(-horizontalDis * (w - 1), horizontalDis * (w - 1), w).tolist()
    zTranslations = np.linspace(-horizontalDis * (l - 1), horizontalDis * (l - 1), l).tolist()
    yTranslations = np.linspace(hStart, hStart + vDis * (h - 1), h).tolist()

    np.random.seed(123456)
    yGap = yTranslations[1] - yTranslations[0]
    for (x, y, z) in product(xTranslations, yTranslations, zTranslations):
        randI = np.random.randint(0,4)
        modelNew = getModel(randI)
        iLayer = int(np.round((y - hStart) / yGap))

        if iLayer >= layerStartToRelaseSequentially:
            modelNew["frameToAppear"] = (iLayer - layerStartToRelaseSequentially) * releaseFrameGap
            y = hStart + yGap * layerStartToRelaseSequentially
        else:
            modelNew["frameToAppear"] = 0

        modelNew['Density'] = 50
        modelNew['translation'][0] = modelNew['translation'][0] + x
        modelNew['translation'][1] = modelNew['translation'][1] + y
        modelNew['translation'][2] = modelNew['translation'][2] + z

        newModels.append(modelNew)



if __name__ == '__main__':
    genFileName = Path(__file__).stem

    machineName = "AnkaPC00"
    # run = False
    run = True

    initialElevation = 5
    elevation = 1.
    r = 2.25

    fixedPoints = []
    models = []

    # make sure you have enough memory, and turn off the renderer
    # numLayers = 72
    # gui = False
    # outputExt = "bin"

    numLayers = 2
    gui = True
    outputExt = "ply"

    addModelsCuboid(models, 10, 12, 12, numLayers, layerStartToRelaseSequentially=100)

    modelsInfo = {
        "Models":models
    }
    parameters = getPhysicsParametersForTest(parametersExample)
    # parameters["PhysicsParams"]["numFrames"] = 300
    parameters["PhysicsParams"]["numFrames"] = 10000
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["worldBounds"] = [
      [  0,  0.0,    0 ],
      [ 0,   20.0,   0  ]
    ]
    recoveryState = None
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True
    parameters["CollisionParams"]["restPoseCloestPoint"] = True

    parameters["PhysicsParams"]["ccdBVHRebuildSteps"] = 100
    parameters["PhysicsParams"]["dcdTetMeshSceneBVHRebuildSteps"] = 100
    parameters["PhysicsParams"]["dcdSurfaceSceneBVHRebuildSteps"] = 100
    parameters["PhysicsParams"]["outputRecoveryStateStep"] = 25
    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.1

    recoveryState = None
    parameters["PhysicsParams"]["outputExt"] = outputExt

    parameters["PhysicsParams"]["useGPU"] = True
    parameters["PhysicsParams"]["collisionStiffness"] = 1e5
    parameters["PhysicsParams"]["numSubsteps"] = 2
    parameters["PhysicsParams"]["iterations"] = 50
    parameters["PhysicsParams"]["binaryModeVisualizationSteps"] = 100

    parameters["PhysicsParams"]["useAccelerator"]=True
    # parameters["PhysicsParams"]["useAccelerator"]=False
    parameters["PhysicsParams"]["acceleratorRho"] = 0.9

    parameters["PhysicsParams"]["saveOutputs"] = False # turn on to output files
    parameters["ViewerParams"] = {
        "enableViewer": gui
    }

    experimentName = "Test_Many_HybridModels_2steps_80Iters_sequential"
    cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=machines[machineName]["binaryFile"], runCommand=run, recoverState=recoveryState)

    experimentName = "Test_Many_HybridModels_2steps_80Iters_sequential_boxDown"
    parameters["PhysicsParams"]["worldBounds"] = [
      [  -18.21156883239746,
        - 10.0, #lower by 10
        -18.280155181884766],
      [ 18.141450881958008,
        409.107666015625,
        18.452653884887695]
    ]
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = False
    cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=machines[machineName]["binaryFile"], runCommand=run, recoverState=recoveryState)


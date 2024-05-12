from M01_Parameters import *
from M02_GenRunningParameters import *
from itertools import product
import numpy as np


def addModelsCuboid(newModels, hStart, w, l, h, vDis = 35, horizontalDis=5.66, rotation=False,
                    rotationAngle = np.pi/3, releaseFrameGap = 75):
    maxVelocityMagnitude = 35
    layerStartToRelaseSequentially = 2
    xTranslations = np.linspace(-horizontalDis * (w - 1), horizontalDis * (w - 1), w).tolist()
    zTranslations = np.linspace(-horizontalDis * (l - 1), horizontalDis * (l - 1), l).tolist()
    yTranslations = np.linspace(hStart, hStart + vDis * (h - 1), h).tolist()
    if h >=2:
        yGap = yTranslations[1] - yTranslations[0]
    else:
        yGap = 1
    np.random.seed(123456)

    for (x, y, z) in product(xTranslations, yTranslations, zTranslations):
        randI = np.random.randint(0,4)
        modelNew = getModelInfoSquishyBallLowRes(modelExample)
        iLayer = int(np.round((y - hStart) / yGap))

        if rotation:

            theta = iLayer * rotationAngle
            R = np.array([
                [np.cos(theta), -np.sin(theta)],
                [np.sin(theta), np.cos(theta)],
                ]
            )

            xz = np.array([[x],[z]])
            xz = R @ xz

            x = float(xz[0])
            z = float(xz[1])

        if iLayer >= layerStartToRelaseSequentially:
            modelNew["frameToAppear"] = (iLayer - layerStartToRelaseSequentially) * releaseFrameGap
            y = hStart + yGap * layerStartToRelaseSequentially
        else:
            modelNew["frameToAppear"] = 0

        modelNew['Density'] = 20
        modelNew['translation'][0] = modelNew['translation'][0] + x
        modelNew['translation'][1] = modelNew['translation'][1] + y
        modelNew['translation'][2] = modelNew['translation'][2] + z

        modelNew['scale'] = [0.3, 0.3, 0.3]
        modelNew['initialVelocity'] = [0, -6, 0]
        modelNew["exponentialVelDamping"] = 0.9
        modelNew["constantVelDamping"] = 0.01
        modelNew['frictionDynamic'] = 0.03
        modelNew["maxVelocityMagnitude"] = maxVelocityMagnitude

        modelNew['dampingHydrostatic'] = 1e-7
        modelNew["dampingDeviatoric"] = 1e-7

        newModels.append(modelNew)

if __name__ == '__main__':
    genFileName = Path(__file__).stem

    machineName = "AnkaPC00"
    # run = False
    run = True

    models = []

    # numLayers = 24
    # gui = False

    numLayers = 1
    gui = True

    addModelsCuboid(models, 40, 3, 3, numLayers, rotation=True)

    teapotContainer = getModelInfoTeapotContainer(modelExample)
    teapotContainer['scale'] = [56, 56, 56]

    models.append(teapotContainer)

    modelsInfo = {
        "Models": models
    }

    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 10000
    # parameters["PhysicsParams"]["numFrames"] = 1
    parameters["PhysicsParams"]["worldBounds"] = [[-10, 0, -10], [10, 20, 10]]

    parameters["PhysicsParams"]["gravity"] = [0.0, -10.0, 0.0]
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True
    parameters["PhysicsParams"]["boundaryFrictionStatic"] = 0.0
    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.0
    parameters["PhysicsParams"]["ccdBVHRebuildSteps"] = 25
    parameters["PhysicsParams"]["dcdBVHRebuildSteps"] = 25
    parameters["PhysicsParams"]["dcdSurfaceSceneBVHRebuildSteps"] = 100

    parameters["PhysicsParams"]["useGPU"] = True
    # parameters["PhysicsParams"]["useGPU"] = False

    parameters["CollisionParams"]["restPoseCloestPoint"] = True
    parameters["PhysicsParams"]["collisionStiffness"] = 5e7

    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["outputRecoveryStateStep"] = 25

    parameters["PhysicsParams"]["useAccelerator"]=True
    # parameters["PhysicsParams"]["useAccelerator"]=False
    parameters["PhysicsParams"]["acceleratorRho"] = 0.92

    parameters["PhysicsParams"]["numSubsteps"] = 5
    parameters["PhysicsParams"]["iterations"] = 40
    parameters["PhysicsParams"]["outputExt"] = "bin"
    parameters["PhysicsParams"]["binaryModeVisualizationSteps"] = 100
    # parameters["PhysicsParams"]["binaryModeVisualizationSteps"] = 50
    parameters["PhysicsParams"]["collisionOffHeight"] = 60

    # parameters["PhysicsParams"]["saveOutputs"] = False # turn on to output files
    parameters["ViewerParams"] = {
        "enableViewer": gui
    }
    experimentName = "Test_216Ball_step5_iter40_largerVDis_squentialRelease"
    recoveryState = r'E:\Data2\VBDSimulation\S05_SquishyBalls_Many_DropToTeapot_sequential\Test_216Ball_step5_iter40_largerVDis_squentialRelease\RecoveryStates\A00000175.json'
    cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=machines[machineName]["binaryFile"],
                                runCommand=run, recoverState=recoveryState)
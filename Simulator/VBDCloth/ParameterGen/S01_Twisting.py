import copy

from M01_Parameters import *
from M02_GenRunningParameters import *
import numpy as np

if __name__ == "__main__":
    genFileName = Path(__file__).stem
    machineName = "PC00"
    # run = False
    run = True

    gui = True
    # binaryFile = r'F:\Code\Projects\Graphics\P15_VBDCloth\Release\P15_VBDCloth.exe'
    binaryFile = machines[machineName]["binaryFile"]

    numberFrames = 1000
    clothSize = 50

    # numberFrames = 15

    contactStiffness = 1e5
    contactDis = 0.2
    queryDis = 0.3
    relaxiation = 0.42
    bendingStiffness = 20
    # experimentName = "Test_SquareCloth"
    model = getModelInfoTestCloth(modelExample, N=clothSize, suffix="_horizontal")

    # fixedPoints = [0, 870]
    # model["fixedPoints"] = fixedPoints
    # model["useParallel"] = True
    model["miu"] = 1e4
    model["lambda"] = 1e4
    model["dampingStVK"] = 3e-5
    model["density"] = 200 * 1e-4
    model["bendingStiffness"] = bendingStiffness

    model["dampingBending"] = 2e-4
    # model["rotation"] = [1.7599884 ,  0, 0]
    model["rotation"] = [0 ,  0, 0]
    model["fixedPoints"] = []

    model["initializationType"] = 2

    modelsInfo = {"Models": [model
                             ]}

    # model["miu"] = 1e5
    # model["lambda"] = 1e5
    parameter = getPhysicsParametersForTest(parametersExample)
    parameter["PhysicsParams"]["timeStep"] = 0.01666666
    parameter["PhysicsParams"]["numSubsteps"] = 10
    parameter["PhysicsParams"]["iterations"] = 5
    parameter["PhysicsParams"]["gravity"] = [0.0, 0.0, 0.0]
    parameter["PhysicsParams"]["useLineSearch"] = False

    parameter["PhysicsParams"]["numFrames"] = numberFrames

    parameter["PhysicsParams"]["outputRecoveryStateStep"] = 10
    parameter["PhysicsParams"]["contactStiffness"] = contactStiffness
    parameter["PhysicsParams"]["contactRadius"] = contactDis
    parameter["PhysicsParams"]["convergenceAvgNormThres"] = 10
    parameter["PhysicsParams"]["conservativeStepRelaxation"] = relaxiation
    parameter["PhysicsParams"]["outputStatistics"] = False

    parameter["PhysicsParams"]["handleCollision"] = True

    parameter["PhysicsParams"]["saveIntermediateResults"] = False

    parameter["ContactDetectorParams"]["maxQueryDis"] = queryDis


    parameter["ViewerParams"] ={
        "enableViewer" : gui
    }

    # parameter["PhysicsParams"]["useParallel"] = False

    leftSide = [clothSize-1 + i*clothSize for i in range(clothSize)]
    rightSide = [0 + i*clothSize for i in range(clothSize)]
    # model["fixedPoints"] = leftSide + rightSide

    deformers = [
        {
            "DeformerName": "Rotator",
            "center": [0, 0, 0],
            "axis": [1, 0, 0],
            "angularVelocity": np.pi / 4,
            "selectedMeshes": [0],
            "selectedVertices": [leftSide],
            "rotationEndTime":60
        },
        {
            "DeformerName": "Rotator",
            "center": [0, 0, 0],
            "axis": [-1, 0, 0],
            "angularVelocity": np.pi / 4,
            "selectedMeshes": [0],
            "selectedVertices": [rightSide],
            "rotationEndTime": 60
        },
    ]
    parameter["Deformers"] = deformers

    experimentName = f"C{clothSize}_Twist_ContactR{contactDis}_CStif{contactStiffness}_bd{bendingStiffness}_damping{model['dampingStVK']}_relaxiation{relaxiation}"

    recoveryState = None

    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameter,
        exeName=binaryFile,
        recoverState=recoveryState,
        runCommand=run,
    )

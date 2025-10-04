import copy

from M01_Parameters import *
from M02_GenRunningParameters import *
import numpy as np
import os

if __name__ == "__main__":
    genFileName = Path(__file__).stem
    machineName = "PC00"
    # run = False
    run = True

    gui = True
    # binaryFile = r'F:\Code\Projects\Graphics\P15_VBDCloth\Release\P15_VBDCloth.exe'
    binaryFile = machines[machineName]["binaryFile"]

    numberFrames = 20000
    clothSize = 100

    # numberFrames = 15

    contactStiffness = 1e6
    contactDis = 0.2
    queryDis = 0.3
    relaxiation = 0.42
    thickness = 1e-4
    bendingStiffness = 50
    pullSpeed = 10
    # experimentName = "Test_SquareCloth"
    model = getModelKnot(modelExample)

    # Get the folder where the script is located
    current_path = os.path.abspath(__file__)
    folder = os.path.dirname(current_path)

    model["path"] = os.path.join(folder, 'KnotGenerator\Data\knot_V2_2400x20_restShape.obj')
    initialShape = r'D:\Code\Graphics\Gaia-Dev\ParameterGen\Python_ParameterGen_VBDCloth\KnotGenerator\Data\knot_V2_2400x20.obj'

    model["verticesColoringCategoriesPath"] = model["path"] + ".vertexColoring.withBendingEnergy.json"

    numVerts = 48000
    numHandles = 20

    # fixedPoints = [0, 870]
    # model["fixedPoints"] = fixedPoints
    # model["useParallel"] = True
    model["miu"] = 1e5
    model["lambda"] = 1e5
    model["dampingStVK"] = 1e-7
    model["density"] = 200 * 1e-4
    model["bendingStiffness"] = bendingStiffness
    model["initialState"] = initialShape
    model["exponentialVelDamping"] = 1
    model["scale"] = [100,100,100]
    model["scaleInitialState"] = True

    model["dampingBending"] = 1e-7
    # model["rotation"] = [1.7599884 ,  0, 0]
    model["rotation"] = [0 ,  0, 0]
    model["fixedPoints"] = []

    model["initializationType"] = 3
    model["use3DRestpose"] = True
    modelsInfo = {"Models": [model
                             ]}

    # model["miu"] = 1e5
    # model["lambda"] = 1e5
    parameter = getPhysicsParametersForTest(parametersExample)
    parameter["PhysicsParams"]["timeStep"] = 0.01666666
    parameter["PhysicsParams"]["numSubsteps"] = 10
    parameter["PhysicsParams"]["iterations"] = 10
    parameter["PhysicsParams"]["gravity"] = [0.0, 0.0, 0.0]
    parameter["PhysicsParams"]["useLineSearch"] = False

    parameter["PhysicsParams"]["numFrames"] = numberFrames

    parameter["PhysicsParams"]["outputRecoveryStateStep"] = 10
    parameter["PhysicsParams"]["outputStatistics"] = False
    parameter["PhysicsParams"]["usePreconditioner"] = True
    parameter["PhysicsParams"]["convergenceAvgNormChangeThres"] = 1e-3
    parameter["PhysicsParams"]["contactStiffness"] = contactStiffness
    parameter["PhysicsParams"]["thickness"] = thickness
    parameter["PhysicsParams"]["contactRadius"] = contactDis
    parameter["PhysicsParams"]["convergenceAvgNormThres"] = 10
    parameter["PhysicsParams"]["conservativeStepRelaxation"] = relaxiation
    parameter["PhysicsParams"]["minStepSizeGD"] = 1e-3
    parameter["PhysicsParams"]["outputStatistics"] = False

    parameter["PhysicsParams"]["handleCollision"] = True

    parameter["PhysicsParams"]["saveIntermediateResults"] = False

    parameter["ContactDetectorParams"]["maxQueryDis"] = queryDis


    parameter["ViewerParams"] ={
        "enableViewer" : gui
    }

    # parameter["PhysicsParams"]["useParallel"] = False

    leftSide = [i for i in range(numHandles)]
    rightSide = [i for i in range(numVerts-numHandles, numVerts)]


    # leftSide = json.load(open(r"W:\02_Graphics\Gaia-Dev\ParameterGen\Python_ParameterGen_VBDCloth\KnotGenerator\Data\knot_leftend.json"))
    # rightSide = json.load(open(r"W:\02_Graphics\Gaia-Dev\ParameterGen\Python_ParameterGen_VBDCloth\KnotGenerator\Data\knot_rightend.json"))
    # model["fixedPoints"] = leftSide + rightSide

    deformers = [
        {
            "DeformerName": "Translator",
            "selectedMeshes": [0],
            "selectedVertices": [leftSide],
            "speed":[-pullSpeed,0,0]
        },
        {
            "DeformerName": "Translator",
            "selectedMeshes": [0],
            "selectedVertices": [rightSide],
            "speed": [pullSpeed, 0, 0]
        },
    ]
    parameter["Deformers"] = deformers

    experimentName = f"Knot_ContactR{contactDis}_CStif{contactStiffness}_bd{bendingStiffness}_damping{model['dampingStVK']}_relaxiation{relaxiation}_pull{pullSpeed}"

    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameter,
        exeName=binaryFile,
        recoverState=None,
        runCommand=run,
    )

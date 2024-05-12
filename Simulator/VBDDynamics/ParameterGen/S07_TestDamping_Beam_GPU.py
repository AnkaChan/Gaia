from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np

if __name__ == "__main__":
    genFileName = Path(__file__).stem

    machineName = "AnkaPC00"
    binaryFile = machines[machineName]["binaryFile"]
    run = True
    model = getModelInfoTest_Beam(modelExample)
    model["fixedPoints"] = [3,4,5,17,18,19,21,46,57]
    modelsInfo = {"Models": [model]}
    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 1000
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["worldBounds"] = [[-100, -100, -100], [100, 100, 100]]
    parameters["PhysicsParams"]["numSubsteps"] = 10
    parameters["PhysicsParams"]["iterations"] = 10
    parameters["PhysicsParams"]["gravity"] = [0.0, -10.0, 0.0]
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True
    parameters["PhysicsParams"]["useDouble3x3"] = False
    parameters["PhysicsParams"]["useGPU"] = True
    parameters["PhysicsParams"]["collisionSolutionType"] = 1

    parameters["PhysicsParams"]["saveOutputs"] = True
    parameters["ViewerParams"] ={
        "enableViewer" : True
    }

    recoveryState = None

    experimentName = "Test_Beam_2e-5_2e-5"
    model["dampingDeviatoric"] = 2e-5
    model["dampingHydrostatic"] = 2e-5
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
        recoverState=recoveryState,
    )

    experimentName = "Test_Beam_5e-5_5e-5"
    model["dampingDeviatoric"] = 5e-5
    model["dampingHydrostatic"] = 5e-5
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
        recoverState=recoveryState,
    )

    experimentName = "Test_Beam_2e-5_0"
    model["dampingDeviatoric"] = 2e-5
    model["dampingHydrostatic"] = 0
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
        recoverState=recoveryState,
    )

    experimentName = "Test_Beam_0_2e-5"
    model["dampingDeviatoric"] = 0
    model["dampingHydrostatic"] = 2e-5
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
        recoverState=recoveryState,
    )


from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np

if __name__ == "__main__":
    genFileName = Path(__file__).stem

    machineName = "AnkaPC00"
    binaryFile = machines[machineName]["binaryFile"]

    run = True
    # run = False

    epsV = 1-3

    cube1 = getModelInfoCubeLowRes(modelExample)
    cube2 = getModelInfoCubeLowRes(modelExample)

    modelsInfo = {"Models": [cube1, cube2]}

    parameters = getPhysicsParametersForTest(parametersExample)

    experimentName = "TestFriction_2Cube_no_friction"
    parameters["PhysicsParams"]["gravity"] = [-5.0, -8.66, 0.0]
    parameters["PhysicsParams"]["numSubsteps"] = 4
    parameters["PhysicsParams"]["iterations"] = 60
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["collisionSolutionType"] = 1
    parameters["PhysicsParams"]["worldBounds"] = [[-100.0, -100.0, -100.0], [100.0, 100.0, 100.0]]
    parameters["PhysicsParams"]["epsV"] = 1e-2
    parameters["PhysicsParams"]["useGPU"] = True
    # parameters["PhysicsParams"]["useGPU"] = False

    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = True
    cube1["translation"] = [0.0, -2.01, 0.0]
    cube1['dampingHydrostatic'] = 1e-7
    cube1["dampingDeviatoric"] = 1e-7
    cube1["frictionDynamic"] = 0.0
    cube1["scale"] = [10.0,1.0,2.0]
    cube1["epsV"] = epsV
    cube1["fixedPoints"] = [i for i in range(100)]

    cube2["translation"] = [5.0, 0.0, 0.0]
    cube2['dampingHydrostatic'] = 1e-7
    cube2["dampingDeviatoric"] = 1e-7
    cube2["dampingGamma"] = 0.0
    cube2["frictionDynamic"] = 0.0
    cube2["epsV"] = epsV
    
    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.0
    parameters["PhysicsParams"]["activeCollisionListPreAllocationRatio"] = 2.0
    parameters["PhysicsParams"]["numFrames"] = 360

    parameters["PhysicsParams"]["saveOutputs"] = True
    parameters["ViewerParams"] ={
        "enableViewer" : True
    }

    parameters["PhysicsParams"]["useAccelerator"]=True
    # parameters["PhysicsParams"]["useAccelerator"]=False
    parameters["PhysicsParams"]["acceleratorRho"] = 0.95


    experimentName = "TestFriction_2Cube_friction_0.1_Iter10"
    cube1["frictionDynamic"] = 0.1
    cube2["frictionDynamic"] = 0.1
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
    )
    experimentName = "TestFriction_2Cube_friction_0.5_Iter10"
    cube1["frictionDynamic"] = 0.5
    cube2["frictionDynamic"] = 0.5
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
    )

    experimentName = "TestFriction_2Cube_friction_0.6"
    cube1["frictionDynamic"] = 0.6
    cube2["frictionDynamic"] = 0.6
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
    )

    experimentName = "TestFriction_2Cube_friction_0.9"
    cube1["frictionDynamic"] = 0.9
    cube2["frictionDynamic"] = 0.9
    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameters,
        exeName=binaryFile,
        runCommand=run,
    )



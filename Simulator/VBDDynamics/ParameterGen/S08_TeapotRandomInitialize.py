from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np


def randomlyInitialize(initialStateFile, outFile):
    # np.random.seed(123154)
    np.random.seed(136)

    data = json.load(open(initialStateFile))
    pts = np.array(data["meshesState"][0]["position"])

    center = np.array([0, 0, 0])
    radius = 0.5

    numPts = pts.shape[0]

    randPts = 2 * np.random.random((numPts, 3)) - np.array([1, 1, 1])

    randPts = randPts / np.linalg.norm(randPts, axis=1)[..., None]

    randPts = randPts * radius + center

    data["meshesState"][0]["position"] = randPts.tolist()

    json.dump(data, open(outFile, "w"))

    return


if __name__ == "__main__":
    genFileName = Path(__file__).stem

    initialState = r"Data\S08_TeapotRandomInitialize\A00000000.json"
    outRandomInitialState = r"Data\S08_TeapotRandomInitialize\A00000000_rand.json"

    randomlyInitialize(initialState, outRandomInitialState)

    machineName = "AnkaPC00"
    # run = False
    binaryFile = machines[machineName]["binaryFile"]
    run = True

    models = []
    numModels = 1

    fixedPoints = [0]
    model = getModelInfoTeapot(modelExample)
    model["miu"] = 2e5
    model["lmbd"] = 1e6
    model["fixedPoints"] = fixedPoints

    model["translation"] = [0, 3, 0]
    model["exponentialVelDamping"] = 0.0
    model["constantVelDamping"] = 0.0
    models.append(model)

    modelsInfo = {"Models": models}

    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 1000
    parameters["PhysicsParams"]["worldBounds"] = [[-10, -10, -10], [10, 20, 10]]

    recoveryState = None

    experimentName = "Test_1_teapot"
    parameters["PhysicsParams"]["gravity"] = [0.0, -10.0, 0.0]
    parameters["CollisionParams"]["allowDCD"] = False
    parameters["CollisionParams"]["allowCCD"] = False
    parameters["PhysicsParams"]["numSubsteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 100

    setToOutputEverySubsteps(parameters)

    parameters["PhysicsParams"]["useGPU"] = True
    # parameters["PhysicsParams"]["useGPU"] = False
    parameters["PhysicsParams"]["outputRecoveryStateStep"] = 10

    parameters["PhysicsParams"]["gravity"] = [0.0, 0.0, 0.0]

    recoveryState = os.path.abspath(outRandomInitialState)
    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.5

    # cmd = genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, exeName=binaryFile,
    #                             runCommand=run, recoverState=recoveryState)
    parameters["PhysicsParams"]["collisionSolutionType"] = 1
    parameters["PhysicsParams"]["numFrames"] = 500
    parameters["PhysicsParams"]["numSubsteps"] = 10
    parameters["PhysicsParams"]["iterations"] = 10
    parameters["PhysicsParams"]["timeStep"] = 1 / 60
    parameters["PhysicsParams"]["stepSizeGD"] = 0
    setToOutputEverySubsteps(parameters)
    parameters["PhysicsParams"]["debugVerboseLvl"] = 1
    model["exponentialVelDamping"] = 1.0
    model["constantVelDamping"] = 0.0

    experimentName = "Test_1_teapot_randomInitialization_damping_0"
    model['dampingHydrostatic'] = 0
    model["dampingDeviatoric"] = 0

    parameters["PhysicsParams"]["saveOutputs"] = True
    parameters["ViewerParams"] ={
        "enableViewer" : True
    }

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


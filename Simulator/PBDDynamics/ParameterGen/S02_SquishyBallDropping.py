import copy
import json
import os, glob
from copy import deepcopy
from itertools import product
import numpy as np
from M01_Parameters import *
from M02_GenRunningParameters import *
from M03_ParametersCreator import *


if __name__ == '__main__':
    genFileName = Path(__file__).stem
    machineName = "MyEnv"
    outputDirId = 0
    run = True

    # numberOfBalls = 16
    numberOfBalls = 4
    vDis = 12.5
    hDis = 7.5
    initialVel = 12.5
    radius = 2.5
    hMin = 20
    hMax = hMin + vDis * numberOfBalls
    numThetaSteps = 4
    numHSteps = int(numberOfBalls / numThetaSteps)

    if numHSteps > 1:
        hStepSize = (hMax - hMin) / (numHSteps - 1)
    else:
        hStepSize = (hMax - hMin)
    hs = np.linspace(hMin, hMax, numHSteps).tolist()
    thetas = np.linspace(0, 2 * np.pi, numThetaSteps, endpoint=False).tolist()

    newModels = []

    for (iH, iTheta) in product(range(numHSteps), range(numThetaSteps)):
        theta = thetas[iTheta]
        h = hs[iH]

        model = getModelInfoSquishyBallLowResHollowMassSpring(modelExample)
        model['translation'] = [radius * np.sin(theta), h + hStepSize * (iTheta / numThetaSteps),
                                   radius * np.cos(theta)]
        model["initialVelocity"] = [0, -initialVel, 0]
        model["density"] = 3
        model["friction_static"] = 0.08
        model["friction_dynamic"] = 0.06
        model["exponentialVelDamping"] = 0.95
        model["constantVelDamping"] = 0.02

        model['devCompliance'] = 5e-4
        model['devDamping'] = 1e-3

        newModels.append(model)

    parameters = getPhysicsParametersoHighResSquishyBall(parametersExample)
    setBVHUpdateStepsForComplicateScene(parameters)
    parameters["PhysicsParams"]["collisionDetectionSubSteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 3
    parameters["PhysicsParams"]["boundary_friction_static"] = 0.0
    parameters["PhysicsParams"]["boundary_friction_dynamic"] = 0.0

    parameters["PhysicsParams"]["gravity"] = [0,0,0]
    parameters["PhysicsParams"]["useBowlGround"] = False
    parameters["PhysicsParams"]["usePlaneGround"] = False
    parameters["PhysicsParams"]["worldBounds"] = [
        [-150, -150, -150],
        [150, 150, 150]
    ]
    parameters["CollisionParams"]["checkFeasibleRegion"] = True
    parameters["CollisionParams"]["rayTriIntersectionEPSILON"] = 1e-10

    parameters["PhysicsParams"]["outputStatistics"] = True
    setBVHUpdateStepsForComplicateScene(parameters)

    modelsInfo = {
        "Models":newModels
    }
    #
    experimentName = f"{numberOfBalls}LowResSquishyBallColliding_SeeNumCollisions"
    genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, recoverState=None,
                          otherArguments=None, outputDirId=outputDirId, runCommand=run)



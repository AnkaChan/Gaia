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


    model = getModelInfo200MSpaghetti(modelExample)

    model["initialVelocity"] = [0, 0, 0]
    model["density"] = 1000
    model["friction_static"] = 0.08
    model["friction_dynamic"] = 0.06
    model["exponentialVelDamping"] = 0.95
    model["constantVelDamping"] = 0.02

    model['devCompliance'] = 5e-4
    model['devDamping'] = 1e-3

    parameters = getPhysicsParametersElasticRod(parametersExample)
    setBVHUpdateStepsForComplicateScene(parameters)
    parameters["PhysicsParams"]["collisionDetectionSubSteps"] = 1
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
        "Models":[model]
    }
    #
    experimentName = "200MNoodle"
    genRunningParameters2(machineName, genFileName, experimentName, modelsInfo, parameters, recoverState=None,
                          otherArguments=None, outputDirId=outputDirId, runCommand=run)



"""
Models Reference:
  "path": "${REPO_ROOT}\\Data\\mesh_models\\ModelsDifferentRes\\sphere-tentacles-v3\\Hollow\\NoPreservingSurface\\sphere2-tentacles1-center85.t",
  "scale": [
    0.3,
    0.3,
    0.3
  ],
  "translation": [
    0.0,
    20.0,
    2.5
  ],
  "Density": 3,
  "DevCompliance": 0.0005,
  "VolCompliance": 0,
  "Damping": 0.001,
  "VelDmping": 0.68,
  "coloredCategoryPath": "${REPO_ROOT}\\Data\\mesh_models\\ModelsDifferentRes\\sphere-tentacles-v3\\Hollow\\NoPreservingSurface\\sphere2-tentacles1-center85.t.coloring.json",
  "edgeColoredCategoryPath": "${REPO_ROOT}\\Data\\mesh_models\\ModelsDifferentRes\\sphere-tentacles-v3\\Hollow\\NoPreservingSurface\\sphere2-tentacles1-center85.t.edgeColoring.json",
  "translationBeforeScaling": [
    0,
    0,
    0
  ],
  "noGravZone": true,
  "initialVelocity": [
    0,
    -12.5,
    0
  ],
  "maxVelocityThreshold": -1,
  "UniformMass": false,
  "friction_static": 0.2,
  "friction_dynamic": 0.20199999999999999,
  "materialName": "MassSpring",
  "noGravZoneThreshold": 60
}

Parameters Reference:
{
  "CollisionParams": {
    "allowCCD": true,
    "allowDCD": true,
    "checkFeasibleRegion": true,
    "checkTetTraverse": true,
    "handleSelfCollision": true,
    "stopTraversingAfterPassingQueryPoint": true,
    "tetrahedralTraverseForNonSelfIntersection": true,
  },
  "General": {
    "checkAndUpdateWorldBounds": true,
    "numTimeSteps": 5000,
    "outputExt": "ply",
    "outputRecoveryState": true,
    "outputRecoveryStateStep": 10,
    "outputStatistics": false,
    "outputT": false,
    "outputVTK": false,
    "saveAllModelsTogether": true,
    "saveSimulationParameters": true,
    "shaderFolderPath": "${REPO_ROOT}/Shader"
  },
  "PhysicsParams": {
    "boundary_friction_dynamic": 0.5,
    "boundary_friction_static": 0.55,
    "bowlCap": false,
    "bowlCenter": [
      0.0,
      25.0,
      0.0
    ],
    "bowlOuterRadius": 25.5,
    "bowlRadius": 25.0,
    "collisionDetectionSubSteps": 1,
    "constantVelDamping": 0.1,
    "debug": false,
    "debugVerboseLvl": 1,
    "doCollDetectionOnlyForFirstIteration": true,
    "doStatistics": false,
    "friction_ground": 0.0,
    "gravity": [
      0.0,
      -10.0,
      0.0
    ],
    "iterations": 2,
    "materialConstraintsInOneStep": false,
    "maxInversionSolveMultiplier": 1.0,
    "maxSurfaceInversionSolveMultiplier": 2.0,
    "numSubsteps": 20,
    "outputIntermediateState": false,
    "perMeshParallelization": false,
    "perTetParallelization": true,
    "restStableDevProjection": false,
    "showSubstepProgress": false,
    "showTimeConsumption": true,
    "smoothSurfaceNormal": true,
    "solveInvertedTets": true,
    "stepInvariantVelDamping": true,
    "timeStep": 0.01666666753590107,
    "useBowlGround": true,
    "usePlaneGround": true,
    "worldBounds": [
      [
        -25.0,
        0.0,
        -25.0
      ],
      [
        25.0,
        50.0,
        25.0
      ]
    ]
  }
}
"""

import numpy as np
from itertools import product
import numpy as np
import os, glob
import json
from copy import deepcopy

def rotationallyStack(modelExample, numModels, modelsStacked, radius, hMin, vdis, numThetaSteps):
    hMax = hMin + numModels * vdis
    numHSteps = int(numModels / numThetaSteps)

    if numHSteps > 1:
        hStepSize = (hMax - hMin) / (numHSteps - 1)
    else:
        hStepSize = (hMax - hMin)
    hs = np.linspace(hMin, hMax, numHSteps).tolist()
    thetas = np.linspace(0, 2 * np.pi, numThetaSteps, endpoint=False).tolist()

    for (iH, iTheta) in product(range(numHSteps), range(numThetaSteps)):
        theta = thetas[iTheta]
        h = hs[iH]

        modelNew = deepcopy(modelExample)
        modelNew['translation'] = [radius * np.sin(theta), h + hStepSize * (iTheta / numThetaSteps),
                                   radius * np.cos(theta)]

        modelsStacked.append(modelNew)

def setSelfCollision(parameters, turnOn = True):
    parameters["CollisionParams"]["handleSelfCollision"] = turnOn


def setCCD(parameters, turnOn):
    parameters["CollisionParams"]["allowCCD"] = turnOn

def setBVHUpdateStepsForComplicateScene(parameters):
    parameters["PhysicsParams"]["dcdTetMeshSceneBVHRebuildSteps"] = 8
    parameters["PhysicsParams"]["dcdSurfaceSceneBVHRebuildSteps"] = 1
    parameters["PhysicsParams"]["ccdBVHRebuildSteps"] = 2

def setBVHUpdateStepsForRealTimeScene(parameters):
    parameters["PhysicsParams"]["dcdTetMeshSceneBVHRebuildSteps"] = 32
    parameters["PhysicsParams"]["dcdSurfaceSceneBVHRebuildSteps"] = 5
    parameters["PhysicsParams"]["ccdBVHRebuildSteps"] = 9
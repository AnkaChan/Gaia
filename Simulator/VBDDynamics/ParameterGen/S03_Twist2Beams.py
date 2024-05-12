from M01_Parameters import *
from M02_GenRunningParameters import *

import numpy as np
import pyvista as pv

if __name__ == "__main__":
    genFileName = Path(__file__).stem
    # machineName = "ZihengLab"
    machineName = "AnkaPC00"
    # machineName = "DataChewer"
    binaryFile = machines[machineName]["binaryFile"]

    run = True
    # run = False
    thinBeam1 = getModelInfoTest_ThinBeamSmoothLowRes(modelExample)
    thinBeam2 = getModelInfoTest_ThinBeamSmoothLowRes(modelExample)

    miu = "5e3"
    lmbd = "1e5"
    collisionStiffness = "2e4"
    thinBeam1["miu"] = float(miu)
    thinBeam1["lmbd"] = float(lmbd)
    thinBeam2["miu"] = float(miu)
    thinBeam2["lmbd"] = float(lmbd)

    thinBeam1["dampingHydrostatic"] = 1e-6
    thinBeam1["dampingDeviatoric"] = 1e-6

    thinBeam2["dampingHydrostatic"] = 1e-6
    thinBeam2["dampingDeviatoric"] = 1e-6

    thinBeam1["frictionDynamic"] = 0.1
    thinBeam2["frictionDynamic"] = 0.1

    thinBeam1["scale"] = [0.1, 0.1, 0.1]
    thinBeam1["translation"] = [0, 0.2, 0]
    thinBeam1["rotation"] = [0, np.pi / 12, 0]
    thinBeam2["miu"] = float(miu)
    thinBeam2["lmbd"] = float(lmbd)
    thinBeam2["scale"] = [0.1, 0.1, 0.1]
    thinBeam2["translation"] = [0, -0.2, 0]
    thinBeam2["rotation"] = [0, -np.pi / 12, 0]

    modelsInfo = {"Models": [thinBeam1, thinBeam2]}

    data = pv.read(
        thinBeam1["vtkPath"].replace("${REPO_ROOT}", machines[machineName]["RepoPath"])
    )
    points = np.array(data.points)
    x_min = np.min(points[:, 0])
    x_max = np.max(points[:, 0])
    x_min_idx = np.where(np.isclose(points[:, 0], x_min, atol=0.5))[0].tolist()
    x_max_idx = np.where(np.isclose(points[:, 0], x_max, atol=0.5))[0].tolist()

    parameters = getPhysicsParametersForTest(parametersExample)
    parameters["PhysicsParams"]["numFrames"] = 100

    parameters["PhysicsParams"]["gravity"] = [0.0, 0.0, 0.0]
    parameters["PhysicsParams"]["numSubsteps"] = 5
    parameters["PhysicsParams"]["iterations"] = 20
    parameters["PhysicsParams"]["checkAndUpdateWorldBounds"] = True
    parameters["PhysicsParams"]["stepSizeGD"] = 0.0
    parameters["PhysicsParams"]["worldBounds"] = [
        [-45, -100.0, -100.0],
        [45, 100.0, 100.0],
    ]
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["allowCCD"] = False
    parameters["PhysicsParams"]["intermediateCollisionIterations"] = 5
    parameters["PhysicsParams"]["checkFeasibleRegion"] = True

    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.0
    parameters["PhysicsParams"]["collisionStiffness"] = float(collisionStiffness)
    parameters["PhysicsParams"]["numFrames"] = 2000
    parameters["PhysicsParams"]["activeCollisionListPreAllocationRatio"] = 2.0
    parameters["PhysicsParams"]["binaryModeVisualizationSteps"] = 10

    deformers = [
        {
            "DeformerName": "Rotator",
            "center": [0, 0, 0],
            "axis": [1, 0, 0],
            "angularVelocity": np.pi / 4,
            "selectedMeshes": [0, 1],
            "selectedVertices": [x_min_idx, x_min_idx],
            "rotationEndTime" :24
        },
        {
            "DeformerName": "Rotator",
            "center": [0, 0, 0],
            "axis": [-1, 0, 0],
            "angularVelocity": np.pi / 4,
            "selectedMeshes": [0, 1],
            "selectedVertices": [x_max_idx, x_max_idx],
            "rotationEndTime": 24
        },
    ]
    recoveryState = None
    experimentName = f"TestThinBeam2Twisting_cpu_{miu}_{lmbd}_{collisionStiffness}_iter{parameters['PhysicsParams']['iterations']}" \
                     f"_damping{thinBeam2['dampingShear']}_{deformers[0]['rotationEndTime']}s"
    parameters["Deformers"] = deformers
    parameters["PhysicsParams"]["useGPU"] = True
    thinBeam1["frictionEpsV"] = 0.3
    thinBeam2["frictionEpsV"] = 0.3
    parameters["PhysicsParams"]["outputExt"] = "bin"

    # turn on to output files
    parameters["PhysicsParams"]["saveOutputs"] = False
    parameters["ViewerParams"] = {
        "enableViewer": True
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
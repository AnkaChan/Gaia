from M01_Parameters import *
from M02_GenRunningParameters import *
import numpy as np

if __name__ == "__main__":
    genFileName = Path(__file__).stem
    machineName = "PC00"
    # run = False
    run = True
    gui = True
    # gui = False
    binaryFile = machines[machineName]["binaryFile"]

    # collider =  os.path.abspath(r"DataGen\Data\SyntheticData\CylinderCollider_tri.obj")
    # colliderParams = copy.deepcopy(modelExample)
    # colliderParams["path"] = collider
    # colliderParams['scale'] = [25,25,25]
    # colliderParams['translation'] = [-18,0,0]
    # colliderParams['rotation'] = [0,0,0]
    # colliderParams["colliderType"] = "TriMeshSequence"
    # colliderParams["meshFiles"] = [
    #     collider
    # ]

    # collider =  os.path.abspath(r"DataGen\Data\SyntheticData\SphereCollider.obj")
    # colliderParams = copy.deepcopy(modelExample)
    # colliderParams["path"] = collider
    # colliderParams['scale'] = [1.5,1.5,1.5]
    # colliderParams['translation'] = [-15,30,0]
    # colliderParams['rotation'] = [0,0,0]
    # colliderParams["colliderType"] = "TriMeshSequence"
    # colliderParams["meshFiles"] = [
    #     collider
    # ]

    numberFrames = 3000

    # numCloth = 50
    numCloth = 5
    disCloth = 0.3

    contactStiffness = 1e5
    maxQueryDis = 0.45
    contactDis = 0.3
    relaxiation = 0.4
    thickness = 1e-4
    bendingStiffness = 10.0
    initType = 2
    res = 70
    accelerationRho = 0.0
    damping = 2e-6

    # experimentName = "Test_SquareCloth"
    modelTemp = getModelInfoTestCloth(modelExample, N=res, suffix="_horizontal")

    randomShift = 0.3

    startHeight = 30
    models = []
    for iCloth in range(numCloth):
        # fixedPoints = [0, 870]
        # model["fixedPoints"] = fixedPoints
        # model["useParallel"] = True
        model = copy.deepcopy(modelTemp)
        model["miu"] = 1e4
        model["lambda"] = 1e4
        model["density"] = 200 * 1e-4
        model["bendingStiffness"] = bendingStiffness
        model["dampingStVK"] = damping
        model["dampingBending"] = damping

        # model["dampingStVK"] = 0
        # model["dampingBending"] = 0
        # model["rotation"] = [1.7599884 ,  0, 0]
        model["rotation"] = [0 ,  0, 0]
        model["translation"] = [0 + randomShift*np.random.randn(), startHeight+ iCloth * disCloth, 0 + randomShift*np.random.randn()]
        # model["fixedPoints"] = [1024, 1124, 1075]
        model["fixedPoints"] = []

        model["initializationType"] = initType

        models.append(model)

    cylinderModel = copy.deepcopy(modelTemp)

    # Get the folder where the script is located
    current_path = os.path.abspath(__file__)
    folder = os.path.dirname(current_path)

    cylinderModel["path"] = os.path.join(folder, r"Data\CylinderCollider_tri.obj")
    cylinderModel["verticesColoringCategoriesPath"] = cylinderModel["path"] + ".vertexColoring.json"
    numVerts = 1224
    cylinderModel["fixedPoints"] = [i for i in range(numVerts)]

    cylinderModel['scale'] = [25, 25, 25]
    cylinderModel['translation'] = [-18,0,0]
    cylinderModel['rotation'] = [0,0,0]
    models.append(cylinderModel)

    modelsInfo = {"Models": models}

    experimentName = f"C50_parallel_tilt{ model['rotation'][0]}_damping{model['dampingStVK']}"
    # model["miu"] = 1e5
    # model["lambda"] = 1e5
    parameter = getPhysicsParametersForTest(parametersExample)

    parameter["PhysicsParams"]["gravity"] = [0.0, -1000.0, 0.0]
    parameter["PhysicsParams"]["worldBounds"] =[
      [
        -10000.0,
        -60.0,
        -10000.0
      ],
      [
        10000.0,
        10000.0,
        10000.0
      ]
    ]

    parameter["PhysicsParams"]["numFrames"] = numberFrames

    parameter["PhysicsParams"]["stepSizeGD"] = 1
    parameter["PhysicsParams"]["outputRecoveryStateStep"] = 10
    parameter["PhysicsParams"]["contactStiffness"] = contactStiffness
    parameter["PhysicsParams"]["contactRadius"] = contactDis
    parameter["ContactDetectorParams"]["maxQueryDis"] = maxQueryDis
    parameter["PhysicsParams"]["thickness"] = thickness
    parameter["PhysicsParams"]["conservativeStepRelaxation"] = relaxiation

    parameter["PhysicsParams"]["usePreconditioner"] = True

    # parameter["PhysicsParams"]["saveIntermediateResults"] = True
    # parameter["PhysicsParams"]["outputStatistics"] = True

    parameter["PhysicsParams"]["evaluateConvergence"] = False
    # parameter["PhysicsParams"]["applyAcceleration"] = True
    parameter["PhysicsParams"]["applyAcceleration"] = False
    parameter["PhysicsParams"]["accelerationRho"] = accelerationRho
    parameter["PhysicsParams"]["outputExt"] = 'bin'
    # parameter["PhysicsParams"]["outputExt"] = 'ply'
    parameter["PhysicsParams"]["binaryModeVisualizationSteps"] = 5
    # parameter["ColliderParams"] = {
    #     "ColliderMeshes":   [colliderParams]
    # }
    parameter["ViewerParams"] ={
        "enableViewer" : gui
    }

    # parameter["PhysicsParams"]["useParallel"] = False

    # leftSide = [49 + i*50 for i in range(50)]
    # rightSide = [0 + i*50 for i in range(50)]
    # model["fixedPoints"] = leftSide + rightSide

    deformers = [

    ]
    parameter["Deformers"] = deformers

    # dt = 0.01
    # steps = 2
    # iters = 20

    dt = 0.016666666666666
    steps = 20
    iters = 15

    # recoveryState = r'E:\Data2\VBD_cloth_Results\S20_LargeScaleTest_collider_ManyCloth\C70_CR0.3_Cs100000.0_bd20.0_step20_iters10_init3_accel0.6_damp2e-06_l100_ball_VBD\RecoveryStates\A00000010.json'
    recoveryState = None
    experimentName = f"C{res}_CR{contactDis}_Cs{contactStiffness}_bd{bendingStiffness}_step{steps}_iters{iters}_init{initType}_accel{accelerationRho}_damp{damping}_l{numCloth}_fix_VBD"
    # experimentName = f"Test_VBD"
    parameter["PhysicsParams"]["timeStep"] = dt
    parameter["PhysicsParams"]["numSubsteps"] = steps
    parameter["PhysicsParams"]["iterations"] = iters
    parameter["PhysicsParams"]["useNewton"] = False
    parameter["PhysicsParams"]["useLineSearch"] = True

    parameter["PhysicsParams"]["handleCollision"] = True
    parameter["PhysicsParams"]["saveAllModelsTogether"] = True

    # cmd = genRunningParameters2(
    #     machineName,
    #     genFileName,
    #     experimentName,
    #     modelsInfo,
    #     parameter,
    #     recoverState=recoveryState,
    #     exeName=binaryFile,
    #     runCommand=run,
    #     commandlineSeparator=' '
    # )

    experimentName = f"C{res}_CR{contactDis}_Cs{contactStiffness}_bd{bendingStiffness}_step{steps}_iters{iters}_init{initType}_accel{accelerationRho}_damp{damping}_l{numCloth}_fix_VBD_run2"
    # experimentName = f"Test_VBD"
    parameter["PhysicsParams"]["timeStep"] = dt
    parameter["PhysicsParams"]["numSubsteps"] = steps
    parameter["PhysicsParams"]["iterations"] = iters
    parameter["PhysicsParams"]["useNewton"] = False
    parameter["PhysicsParams"]["useLineSearch"] = True

    parameter["PhysicsParams"]["handleCollision"] = True
    parameter["PhysicsParams"]["saveAllModelsTogether"] = True

    cmd = genRunningParameters2(
        machineName,
        genFileName,
        experimentName,
        modelsInfo,
        parameter,
        recoverState=recoveryState,
        exeName=binaryFile,
        runCommand=run,
        commandlineSeparator=' '
    )


    dt = 0.01
    steps = 1
    iters = 30
    experimentName = f"C{res}_ContactR{contactDis}_ContactStiffness{contactStiffness}_bd{bendingStiffness}_step{steps}_iters{iters}_init{initType}_newton"
    parameter["PhysicsParams"]["timeStep"] = dt
    parameter["PhysicsParams"]["numSubsteps"] = steps
    parameter["PhysicsParams"]["iterations"] = iters
    parameter["PhysicsParams"]["contactStiffness"] = contactStiffness
    parameter["PhysicsParams"]["contactRadius"] = contactDis
    parameter["ContactDetectorParams"]["maxQueryDis"] = maxQueryDis
    parameter["PhysicsParams"]["thickness"] = thickness
    parameter["PhysicsParams"]["conservativeStepRelaxation"] = relaxiation

    parameter["PhysicsParams"]["useNewton"] = True
    parameter["PhysicsParams"]["useLineSearch"] = True
    parameter["PhysicsParams"]["psd"] = True
    parameter["PhysicsParams"]["debugVerboseLvl"] = 3
    parameter["PhysicsParams"]["iterations"] = iters
    parameter["PhysicsParams"]["handleCollision"] = True
    # parameter["PhysicsParams"]["stepSize"] = 0.5

    # cmd = genRunningParameters2(
    #     machineName,
    #     genFileName,
    #     experimentName,
    #     modelsInfo,
    #     parameter,
    #     exeName=binaryFile,
    #     runCommand=run,
    # )
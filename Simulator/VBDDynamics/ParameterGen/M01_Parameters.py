import copy
import json
from copy import deepcopy
import pathlib

modelExample = json.load(open(str(pathlib.Path(__file__).parent) + "/Examples/Models.json"))["Models"][0]
parametersExample = json.load(open(str(pathlib.Path(__file__).parent) + "/Examples/Parameters.json"))


def setToOutputEverySubsteps(physicsParameters):
    physicsParameters["PhysicsParams"]["timeStep"] = physicsParameters["PhysicsParams"]["timeStep"]  / physicsParameters["PhysicsParams"]["numSubsteps"]
    physicsParameters["PhysicsParams"]["numSubsteps"] = 1

def getModelInfoTestCloth_unitTest2Tets(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}/CPP_Simulation\P10_VBDDynamics\TestData/TwoTets.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e4
    model['lmbd'] = 1e4

    model['density'] = 1000
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_Octupus(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/Octopus3_Inflated/Octopus_Inflated_5k_noIntersection.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_Beam(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/UnitTest/Beam/Beam.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_Beam2(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/UnitTest/Beam2/Beam.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_Beam_1k(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/UnitTest/Beam/Beam_1k.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_ThinBeamHighRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/Subdivided/thinbeam_Subdivided.t"
    model['vtkPath'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/Subdivided/thinbeam_Subdivided.1.vtk"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model


def getModelInfoTest_ThinBeamHighRes2(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/UnitTest/ThinBeam2/ThinBeam2HighRes.t"
    model['vtkPath'] = "${REPO_ROOT}/Data/mesh_models/UnitTest/ThinBeam2/ThinBeam2HighRes.1.vtk"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model


def getModelInfoTest_ThinBeamLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/ThinbeamCentered.t"
    model['vtkPath'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/ThinbeamCentered.1.vtk"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_ThinBeamSmoothLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/ThinbeamCentered_smooth.t"
    model['vtkPath'] = "${REPO_ROOT}/Data/mesh_models/ModelsDifferentRes/thinbeam/ThinbeamCentered_smooth.1.vtk"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTest_BunnyLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = "${REPO_ROOT}/Data/mesh_models/t/bunny_small.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['miu'] = 1e5
    model['lmbd'] = 1e6

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
    else:
        print("Unrecognized material name! Use default material: NeoHookean")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoHighArmadilo(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\t\armadillo.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e6
        model['lmbd'] = 1e7

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "MassSpring"

    return model

def getModelInfoArmadilo15K(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\t\Armadilo_15K.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e6
        model['lmbd'] = 1e7

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "MassSpring"

    return model


def getModelInfoHighArmadiloLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\t\Armadilo_lowres.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e5
        model['lmbd'] = 1e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "MassSpring"

    return model



def getModelInfoWreckingBall(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\MassRatio\WreckingBall.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 1000
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e5
        model['lmbd'] = 1e5

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "MassSpring"

    return model

def getModelInfoWreckingBallLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\MassRatio\WreckingBall_lowRes_bigger.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\MassRatio\WreckingBall_bigger_32Segs.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 1000
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e5
        model['lmbd'] = 1e5

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "MassSpring"

    return model

def getModelInfoCubeLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Cube\Cube.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e5
        model['lmbd'] = 2e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoCubeRandomizedTessellation(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\CubeRandom\Cube.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e5
        model['lmbd'] = 2e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoSquishyBallLowRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\sphere-tentacles-v3\Hollow\sphere2-tentacles1-center85.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 50
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e7
        model['lmbd'] = 1e8

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoSquishyBallHighRes(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\..\Data\TetMesh\Hollow\sphere2-tentacles2_tri_hollow_9.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 50
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e7
        model['lmbd'] = 1e8

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTeapot(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\Teapot\teapot.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 50
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e6
        model['lmbd'] = 2e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoTeapotContainer(modelExample, materialType="NeoHookean", fix=True):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\Containers\TeapotContainer\teapotContainer.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 5000
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e6
        model['lmbd'] = 2e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    if fix:
        numPts = 1165
        model["fixedPoints"] = [i for i in range(numPts)]

    return model


def getModelInfoBox(modelExample, materialType="NeoHookean", fix=True):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\Containers\GlassBox.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 5000
    model['translationBeforeScaling'] = [0, 0.2, 0]

    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 2e7
        model['lmbd'] = 2e8

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    if fix:
        numPts = 16
        model["fixedPoints"] = [i for i in range(numPts)]

    return model

# a sphere centered at (0,0,0) with radius 1
def getModelInfoSphere(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\t\sphere.t"
    model['verticesColoringCategoriesPath'] = model['path'] + ".vertexColoring.json"

    model['density'] = 100
    model['translationBeforeScaling'] = [0, 0, 0]

    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['miu'] = 1e5
        model['lmbd'] = 1e6

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"
    #
    #     model['springStiffness'] = 1e4
    #     model['damping'] = 1e-2

    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getPhysicsParametersForTest(parameterExample):
    parameters = deepcopy(parameterExample)

    parameters["PhysicsParams"]["worldBounds"] = [
        [-10.0, 0.0, -10.0],
        [10.0, 20.0, 10.0]
    ]

    parameters["PhysicsParams"]["gravity"] = [0.0, -10.0, 0.0]

    parameters["PhysicsParams"]["timeStep"] = 0.01666666
    parameters["PhysicsParams"]["numSubsteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 20

    return  parameters



def changeToOutputPerSubStep(parameters):
    parameters["PhysicsParams"]["timeStep"] = parameters["PhysicsParams"]["timeStep"] / parameters["PhysicsParams"]["numSubsteps"]
    parameters["PhysicsParams"]["numSubsteps"] = 1
    return parameters
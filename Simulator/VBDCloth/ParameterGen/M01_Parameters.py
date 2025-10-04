import copy
import json
from copy import deepcopy

modelExample = json.load(open("./Examples/Models.json"))["Models"][0]
parametersExample = json.load(open("./Examples/Parameters.json"))


def setToOutputEverySubsteps(physicsParameters):
    physicsParameters["PhysicsParams"]["timeStep"] = (
        physicsParameters["PhysicsParams"]["timeStep"]
        / physicsParameters["PhysicsParams"]["numSubsteps"]
    )
    physicsParameters["PhysicsParams"]["numSubsteps"] = 1


def getModelInfoTestCloth(
    modelExample, N=30, suffix="", materialType="StVK_triMesh", colorWithBending=True
):
    model = copy.deepcopy(modelExample)
    model["path"] = (
        rf"${{REPO_ROOT}}\Simulator\VBDCloth\ParameterGen\Data\C{N}{suffix}.obj"
    )
    if colorWithBending:
        model["verticesColoringCategoriesPath"] = (
            model["path"] + ".vertexColoring.withBendingEnergy.json"
        )
    else:
        model["verticesColoringCategoriesPath"] = model["path"] + ".vertexColoring.json"

    model["miu"] = 1e4
    model["lambda"] = 1e4

    model["density"] = 200 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"

    return model

def getModelKnot(modelExample, materialType="StVK_triMesh"):
    model = copy.deepcopy(modelExample)
    model["path"] = (
        rf"${{REPO_ROOT}}\ParameterGen\Python_ParameterGen_VBDCloth\KnotGenerator\Data\knot_restShape.obj"
    )
    model["verticesColoringCategoriesPath"] = (
            model["path"] + ".vertexColoring.withBendingEnergy.json"
        )
    model["miu"] = 1e4
    model["lambda"] = 1e4

    model["density"] = 200 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"
    return model



def getModelInfoSuitSquare(
    modelExample, dir, suffix="", materialType="StVK_triMesh", colorWithBending=True
):
    model = copy.deepcopy(modelExample)
    model["path"] = rf"{dir}\square{suffix}.obj"
    if colorWithBending:
        model["verticesColoringCategoriesPath"] = (
            model["path"] + ".vertexColoring.withBendingEnergy.json"
        )
    else:
        model["verticesColoringCategoriesPath"] = model["path"] + ".vertexColoring.json"

    model["miu"] = 1e4
    model["lambda"] = 1e4

    model["density"] = 200 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"

    return model

def getModelInfo_dress_knife(dataFolder,
    modelExample=modelExample, materialType="StVK_triMesh"
):
    model = copy.deepcopy(modelExample)
    model["path"] = (
        rf"{dataFolder}\rest_shape.obj"
    )
    model["verticesColoringCategoriesPath"] = (
        model["path"] + ".vertexColoring.withBendingEnergy.json"
    )
    model["miu"] = 1e5
    model["lambda"] = 1e5

    model["density"] = 174 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"
    return model


def getModelInfoUnitTest1Triangle(modelExample, suffix="", materialType="StVK_triMesh"):
    model = copy.deepcopy(modelExample)
    model["path"] = (
        rf"${{REPO_ROOT}}\ParameterGen\Python_ParameterGen_VBDCloth\Data\SyntheticData\Triangle_1_on_C5.obj"
    )
    model["verticesColoringCategoriesPath"] = model["path"] + ".vertexColoring.json"

    model["miu"] = 1e4
    model["lambda"] = 1e4

    model["density"] = 200 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"

    return model


def getModelInfoTestCloth_unitTest(modelExample, materialType="StVK_triMesh"):
    model = copy.deepcopy(modelExample)
    model["path"] = (
        r"${REPO_ROOT}\/ParameterGen/Python_ParameterGen_APAPCloth/DataGen/Data/SyntheticData/C4.obj"
    )
    model["verticesColoringCategoriesPath"] = model["path"] + ".vertexColoring.json"

    model["miu"] = 1e4
    model["lambda"] = 1e4

    model["density"] = 200 * 1e-4
    if materialType == "StVK_triMesh":
        model["materialName"] = "StVK_triMesh"
        model["bendingStiffness"] = 25
    else:
        print("Unrecognized material name!")
        model["materialName"] = "StVK_triMesh"

    return model


def getPhysicsParametersForTest(parameterExample):
    parameters = deepcopy(parameterExample)

    parameters["PhysicsParams"]["gravity"] = [0.0, -1000.0, 0.0]

    parameters["PhysicsParams"]["timeStep"] = 0.01666666
    parameters["PhysicsParams"]["numSubsteps"] = 20
    parameters["PhysicsParams"]["iterations"] = 2

    return parameters

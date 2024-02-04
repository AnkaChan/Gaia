import copy
import json
from copy import deepcopy

modelExample = json.load(open("./Examples/Models.json"))["Models"][0]
parametersExample = json.load(open("./Examples/Parameters.json"))

def getModelInfoOctopusInflatted(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Octopus3_Inflated\Octopus_Inflated_5k_noIntersection.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"

    model['density'] = 50
    model['scale'] = [1, 1, 1]
    model["hasNoGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["exponentialVelDamping"] = 0.75
    model["constantVelDamping"] = 0.0
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.02

    elif materialType == "MassSpring":
        model["materialName"] = "MassSpring"
        model['edgesColoringCategoriesPath'] = model['path'] + ".edgecoloring.json"

        model['springCompliance'] = 1e-4
        model['springDamping'] = 0.02
    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model


def getPhysicsParametersoOctopus(parameterExample):
    parameters = deepcopy(parameterExample)

    parameters["PhysicsParams"]["worldBounds"] = [
        [-10.0, 0.0, -10.0],
        [ 10.0, 20.0, 10.0]
    ]
    parameters["PhysicsParams"]["numSubsteps"] = 20
    parameters["PhysicsParams"]["iterations"] = 2
    return  parameters

def getModelInfoLowResSquishyBall(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\sphere-tentacles-v3\Hollow\sphere2-tentacles1-center85.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 5
    model['scale'] = [0.3, 0.3, 0.3]
    model["hasNoGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["translation"] = [0, 10, 0]
    model["maxVelocityThreshold"] = -1
    model["exponentialVelDamping"] = 0.75
    model["constantVelDamping"] = 0.0
    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.001

    elif materialType == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 5e-4
        model['springDamping'] = 0.001
    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getModelInfoHighResSquishyBall(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\..\Data\TetMesh\Hollow\sphere2-tentacles2_tri_hollow_9.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 5
    model['scale'] = [0.3, 0.3, 0.3]
    model["hasNoGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["translation"] = [0, 10, 0]
    model["maxVelocityThreshold"] = -1
    model["exponentialVelDamping"] = 0.75
    model["constantVelDamping"] = 0.01
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.001

    elif materialType == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 5e-4
        model['springDamping'] = 0.001
    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    return model

def getPhysicsParametersoHighResSquishyBall(parameterExample):
    parameters = deepcopy(parameterExample)

    parameters["PhysicsParams"]["worldBounds"] = [
        [-50.0, 0.0, -50.0],
        [ 50.0, 200.0, 50.0]
    ]
    parameters["PhysicsParams"]["numSubsteps"] = 50
    parameters["PhysicsParams"]["iterations"] = 2

    parameters["CollisionParams"]["feasibleRegionEpsilon"] = 1e-4
    parameters["CollisionParams"]["maxSearchDistanceMultiplier"] = 1.5

    return  parameters

def getModelInfoLowResThinbeam(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)

    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\thinbeam.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\ThinbeamCentered.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 1000
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [0.1, 0.1, 0.1]

    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False

    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05

    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 1e-1
        model['springDamping'] = 0.1
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.01

    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    return model

def getModelInfoLowResCuBoid(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)

    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\thinbeam.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Cuboid\Cuboid.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 1000
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [1, 1, 1]

    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False

    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05

    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 1e-1
        model['springDamping'] = 0.1
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.01

    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    return model

def getModelInfoOverhandKnotFractal(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)

    # modelExample['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\knots\OverhandKnot.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\knots\OverhandKnot_fractal.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['density'] = 2500
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [1, 1, 1]
    model['devCompliance'] = 1e-4
    model['damping'] = 1e-3
    model['volCompliance'] = 0
    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False
    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['devCompliance'] = 1e-4
        model['volCompliance'] = 0
        model['Damping'] = 1e-3
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['volCompliance'] = 0
        model['damping'] = 1e-3

    model["friction_static"] = 0.04
    model["friction_dynamic"] = 0.03

    return model

def getModelInfoHighResCuBoid(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)

    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\thinbeam.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Cuboid\Cuboid_subdivided.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 1000
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [1, 1, 1]

    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False

    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05

    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 1e-1
        model['springDamping'] = 0.1
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.01

    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    return model


def getModelInfoCube(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)
    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\thinbeam.t"
    # model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Cube\Cube_subdivided.t"
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\Cube\Cube.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".tetColoring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 1000
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [1, 1, 1]

    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False

    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05

    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 1e-1
        model['springDamping'] = 0.1
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.01

    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    return model

def getPhysicsParametersOverhandKnot(parameters, parallelizationType='PerTet', boundary = "plane", bowlRadius=30):
    parameters = copy.deepcopy(parameters)

    parameters["CollisionParams"]["allowCCD"] = True
    parameters["CollisionParams"]["allowDCD"] = True
    parameters["CollisionParams"]["stopTraversingAfterPassingQueryPoint"] = False
    parameters["CollisionParams"]["checkFeasibleRegion"] = False

    parameters["PhysicsParams"]["numTimeSteps"] = 5000
    parameters["PhysicsParams"]["outputVTK"] = False

    # step & iterations
    parameters["PhysicsParams"]["numSubsteps"] = 50
    parameters["PhysicsParams"]["collisionDetectionSubSteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 3
    parameters["PhysicsParams"]["doCollDetectionOnlyForFirstIteration"] = True

    # other solver settings
    parameters["PhysicsParams"]["stepInvariantVelDamping"] = True
    parameters["PhysicsParams"]["smoothSurfaceNormal"] = True

    # boundaries
    parameters["PhysicsParams"]["useBowlGround"] = False
    parameters["PhysicsParams"]["usePlaneGround"] = True
    parameters["PhysicsParams"]["worldBounds"] = [
        [-1000, -1000, -1000],
        [1000, 1000, 1000]
    ]

    # Material setting
    parameters["PhysicsParams"]["restStableDevProjection"] = False
    parameters["PhysicsParams"]["materialConstraintsInOneStep"] = False
    parameters["PhysicsParams"]["solveInvertedTets"] = True

    # Friction setting
    parameters["PhysicsParams"]["friction_ground"] = 1000
    parameters["PhysicsParams"]["boundary_friction_static"] = 0.55
    parameters["PhysicsParams"]["boundary_friction_dynamic"] = 0.5

    # Debug Configuration
    parameters["PhysicsParams"]["debug"] = False
    parameters["PhysicsParams"]["showSubstepProgress"] = False
    parameters["PhysicsParams"]["debugLevel"] = 1
    # parameters["PhysicsParams"]["doStatistics"] = True

    parameters["PhysicsParams"]["gravity"] = [0, 0, 0]

    return parameters

def getModelInfoThinbeamSubdivided(modelExample, materialName="MassSpring"):
    model = copy.deepcopy(modelExample)

    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\thinbeam\Subdivided\thinbeam_Subdivided.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoredCategoryPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 25
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 0, 0]
    model['scale'] = [0.1, 0.1, 0.1]

    model["noGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["maxVelocityThreshold"] = -1
    model["uniformMass"] = False

    model["exponentialVelDamping"] = 0.4
    model["constantVelDamping"] = 0.05

    if materialName == "MassSpring":
        model["materialName"] = "MassSpring"
        model['springCompliance'] = 1e-1
        model['springDamping'] = 0.1
    else:
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-3
        model['devDamping'] = 0.01

    model["friction_static"] = 0.16
    model["friction_dynamic"] = 0.12
    return model



def getPhysicsParametersThinBeam(parameters, parallelizationType='PerTet', boundary = "plane", bowlRadius=30):
    parameters = copy.deepcopy(parameters)

    parameters["CollisionParams"]["allowCCD"] = True
    parameters["CollisionParams"]["allowDCD"] = True

    parameters["CollisionParams"]["stopTraversingAfterPassingQueryPoint"] = False
    parameters["CollisionParams"]["checkFeasibleRegion"] = True

    parameters["PhysicsParams"]["numTimeSteps"] = 5000
    parameters["PhysicsParams"]["outputVTK"] = False

    # step & iterations
    parameters["PhysicsParams"]["numSubsteps"] = 100
    parameters["PhysicsParams"]["collisionDetectionSubSteps"] = 1
    parameters["PhysicsParams"]["iterations"] = 2
    parameters["PhysicsParams"]["doCollDetectionOnlyForFirstIteration"] = True
    parameters["PhysicsParams"]['timeStep'] = 0.0333333

    # other solver settings
    parameters["PhysicsParams"]["stepInvariantVelDamping"] = True
    parameters["PhysicsParams"]["smoothSurfaceNormal"] = True

    # boundaries
    parameters["PhysicsParams"]["usePlaneGround"] = True
    parameters["PhysicsParams"]["useBowlGround"] = False

    parameters["PhysicsParams"]["worldBounds"] = [
        [-5.1, -200, -200],
        [5.1, 200, 200]
    ]

    # Parallelization
    if parallelizationType == 'PerTet':
        parameters["PhysicsParams"]["perMeshParallelization"] = False
        parameters["PhysicsParams"]["perTetParallelization"] = True
    elif parallelizationType == 'PerMesh':
        parameters["PhysicsParams"]["perMeshParallelization"] = True
        parameters["PhysicsParams"]["perTetParallelization"] = False
    elif parallelizationType == 'None':
        parameters["PhysicsParams"]["perMeshParallelization"] = False
        parameters["PhysicsParams"]["perTetParallelization"] = False
    else:
        assert False

    # Material setting
    parameters["PhysicsParams"]["restStableDevProjection"] = False
    parameters["PhysicsParams"]["materialConstraintsInOneStep"] = False
    parameters["PhysicsParams"]["solveInvertedTets"] = True

    # Friction setting
    parameters["PhysicsParams"]["friction_ground"] = 1000
    parameters["PhysicsParams"]["boundaryFrictionDynamic"] = 0.12
    parameters["PhysicsParams"]["boundaryFrictionStatic"] = 0.15

    # Debug Configuration
    parameters["PhysicsParams"]["debug"] = False
    parameters["PhysicsParams"]["showSubstepProgress"] = False
    parameters["PhysicsParams"]["debugLevel"] = 1
    # parameters["PhysicsParams"]["doStatistics"] = True

    parameters["PhysicsParams"]["gravity"] = [0, -0.25, 0]



    return parameters

def getModelInfoElasticRod(modelExample, materialType="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\PropellerRubberBand\RubberBand.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"

    model['density'] = 2500
    model["hasNoGravZone"] = False
    model["initialVelocity"] = [0, 0, 0]
    model["translation"] = [0, 10, 0]
    model["maxVelocityThreshold"] = -1
    model["exponentialVelDamping"] = 0.75
    model["constantVelDamping"] = 0.0
    if materialType == "NeoHookean":
        model["materialName"] = "NeoHookean"
        model['devCompliance'] = 1e-4
        model['devDamping'] = 0.001

    # elif materialType == "MassSpring":
    #     model["materialName"] = "MassSpring"
    #     model['springCompliance'] = 5e-4
    #     model['springDamping'] = 0.001
    else:
        print("Unrecognized material name!")
        model["materialName"] = "NeoHookean"

    model["friction_static"] = 0.12
    model["friction_dynamic"] = 0.10
    return model


def getPhysicsParametersElasticRod(parameterExample):
    parameters = deepcopy(parameterExample)

    parameters["PhysicsParams"]["worldBounds"] = [
        [-50.0, 0.0, -50.0],
        [ 50.0, 200.0, 50.0]
    ]
    parameters["PhysicsParams"]["numSubsteps"] = 60
    parameters["PhysicsParams"]["iterations"] = 6

    parameters["PhysicsParams"]["worldBounds"] = [
        [-1000, -1000, -1000],
        [1000, 1000, 1000]
    ]
    parameters["PhysicsParams"]["numFrames"] = 20000


    return  parameters


def getModelInfoSquishyBallLowResHollowMassSpring(modelExample, modelPreverseSurface=False):
    model = copy.deepcopy(modelExample)

    if modelPreverseSurface:
        model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\sphere-tentacles-v3\Hollow\sphere2-tentacles1-center85.t"

    else:
        model['path'] = r"${REPO_ROOT}\Data\mesh_models\ModelsDifferentRes\sphere-tentacles-v3\Hollow\NoPreservingSurface\sphere2-tentacles1-center85.t"

    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"
    model['density'] = 2
    model['translationBeforeScaling'] = [0, 0, 0]
    model['translation'] = [0, 10, 0]
    model['scale'] = [0.3, 0.3, 0.3]
    model['devCompliance'] = 5e-4
    model['damping'] = 1e-3

    model["noGravZone"] = False
    model["initialVelocity"] = [0, -3, 0]
    model["velDmping"] = 0.68
    model["maxVelocityThreshold"] = -1
    model["UniformMass"] = False

    model["friction_static"] = 0.2
    model["friction_dynamic"] = 0.78 * model["friction_static"]

    model["materialName"] = "MassSpring"

    return  model

def getModelInfo200MSpaghetti(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)

    model['path'] = r"${REPO_ROOT}\..\Data\TetMesh\Spaghetti_1m_lowRes_tri_Stack_201.t"
    model['tetsColoringCategoriesPath'] = model['path'] + ".coloring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"
    model["initialVelocity"] = [ 0, -3, 0]
    model["maxVelocityMagnitude"] = 4.5
    return model

def getModelInfoWreckingBallLowRes(modelExample, materialName="NeoHookean"):
    model = copy.deepcopy(modelExample)
    model["materialName"] = "NeoHookean"
    model['devCompliance'] = 1e-8
    model['volCompliance'] = 1e-8
    model['path'] = r"${REPO_ROOT}\Data\mesh_models\UnitTest\MassRatio\WreckingBall_lowRes_bigger.t"

    model['density'] = 1000

    model['tetsColoringCategoriesPath'] = model['path'] + ".tetColoring.json"
    model['edgesColoringCategoriesPath'] = model['path'] + ".edgeColoring.json"
    # model["initialVelocity"] = [ 0, -3, 0]
    # model["maxVelocityMagnitude"] = 4.5
    return model

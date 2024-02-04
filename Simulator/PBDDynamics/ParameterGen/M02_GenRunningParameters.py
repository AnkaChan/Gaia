import json, subprocess, time
import os
from os.path import join
from pathlib import Path

machines = {
    "MyEnv":{
        "RepoPath":None,
        "ExecutablePath":None,
        "OutputDirs": [None],
    },
    "TestEnv": {
        "RepoPath": r"F:\Code\02_Graphics\Gaia",
        "ExecutablePath": r"F:\Code\Projects\Graphics\PBDDynamics_GAIA\Release\PBDDynamics.exe",
        "OutputDirs": [r"E:\Data2\E_PBDCollision2Results"],
    },
}

def genRunningParameters(outputCmdFile, machineName, experimentPath, modelsFile, parametersFile, experimentName, outputDirId=0,
                         exeName = r"P02_PBDDynamics.exe", recoverState = None, relativeRecoverStatePath=True, otherArguments = None):
    machineInfo = machines[machineName]
    outputDir = machineInfo["OutputDirs"][outputDirId]

    cmd = [
        exeName,
        join(machineInfo["RepoPath"], "ExperimentParameter", experimentPath, modelsFile),
        join(machineInfo["RepoPath"], "ExperimentParameter", experimentPath, parametersFile),
        join(outputDir, experimentPath, experimentName),
        '-R', machineInfo["RepoPath"]
    ]

    if recoverState is not None:
        if relativeRecoverStatePath:
            recoverState = join(outputDir, recoverState)
        cmd.append("-r " + recoverState)

    if "temp" in exeName:
        print("Warning!!! temp exe is used!")

    if otherArguments is not None:
        cmd.extend(otherArguments)
    elif exeName == r"P12_PBDPhysicsWithRotator.exe":
            print(exeName, " requires extra arguments!!!")


    f = open(outputCmdFile, 'w')
    f.write(" ".join(cmd))
    f.write("\ncmd /k")

def genRunningParameters2(machineName, experimentPath, experimentName, modelsInfo, parametersInfo, modelsFile="Models.json", parametersFile="Parameters.json", outputDirId=0,
                          recoverState = None, relativeRecoverStatePath=True, otherArguments = None, runCommand=False, absoluteBinaryPath=False):
    machineInfo = machines[machineName]
    outputDir = machineInfo["OutputDirs"][outputDirId]

    if outputDir is None:
        print("Please write the preferred output path to \Gaia\Simulator\PBDDynamics\ParameterGen\M02_GenRunningParameters.py: machines['MyEnv']['OutputDirs']")
        exit(-1)

    parameterOutputDir = join(".", "Parameters", experimentPath, experimentName)
    os.makedirs(parameterOutputDir, exist_ok=True)

    json.dump(modelsInfo, open(join(parameterOutputDir, modelsFile), 'w'), indent=2)
    json.dump(parametersInfo, open(join(parameterOutputDir, parametersFile), 'w'), indent=2)


    if machineInfo["RepoPath"] is None:
        print("Please write the absolute repo path to \Gaia\Simulator\PBDDynamics\ParameterGen\M02_GenRunningParameters.py: machines['MyEnv']['RepoPath']")
        exit(-1)

    modelsFileRemote = join(machineInfo["RepoPath"], "Simulator", "PBDDynamics", "ParameterGen", "Parameters", experimentPath, experimentName, modelsFile)
    parametersFileRemote = join(machineInfo["RepoPath"], "Simulator", "PBDDynamics", "ParameterGen", "Parameters", experimentPath, experimentName, parametersFile)

    if machineInfo["ExecutablePath"] is None:
        print(
            "Please write the absolute path of the executable to \Gaia\Simulator\PBDDynamics\ParameterGen\M02_GenRunningParameters.py: machines['MyEnv']['ExecutablePath']")
        exit(-1)

    cmd = [
        machineInfo["ExecutablePath"],
        modelsFileRemote,
        parametersFileRemote,
        join(outputDir, experimentPath, experimentName),
        '-R', machineInfo["RepoPath"]
    ]

    if recoverState is not None:
        if relativeRecoverStatePath:
            recoverState = join(outputDir, recoverState)
        cmd.append("-r " + recoverState)

    if otherArguments is not None:
        if isinstance(otherArguments, list):
            cmd.extend(otherArguments)
        elif isinstance(otherArguments, str):
            cmd.append(otherArguments)

    outputCmdFile = join(parameterOutputDir, "Run_" + machineName + ".cmd")
    f = open(outputCmdFile, 'w')
    f.write(" ".join(cmd))
    print("Command: ", "\n".join(cmd))

    f.write("\ncmd /k")

    if runCommand:
        startT = time.time()
        subprocess.run(cmd)
        T = time.time() - startT
        print("Running command took: ", T)


if __name__ == '__main__':
    # generate thinbeam rotation running parameters for DataChewer
    # rotatorArguments = ["-a", "3.14", "-e", "20"]
    # recoveryState = r"Demo23_Twist\ThinBeam_subdivided_NeoHookean\run1\RecoveryStates\A00000590.json"
    # genRunningParameters("../../Bin/Run_DataChewer_ThinBeam_Neohookean.cmd", "DataChewer", "Demo23_Twist\ThinBeam_subdivided_NeoHookean", "Models.json",
    #                      "Parameters.json", "run2_withStop",
    #                      exeName = "P12_PBDPhysicsWithRotator.exe", otherArguments=rotatorArguments, recoverState=recoveryState )

    # generate thinbeam rotation running parameters for DataChewer
    rotatorArguments = ["-a", "3.14", "-e", "20"]
    recoveryState = r"Demo23_Twist\ThinBeam_subdivided_NeoHookean\run1\RecoveryStates\A00000590.json"
    genRunningParameters("../../Bin/Run_AnkaPC01_ThinBeam_Neohookean.cmd", "AnkaPC01", "Demo23_Twist\ThinBeam_subdivided_NeoHookean", "Models.json",
                         "Parameters.json", "run6_withStop_CCD_VelDamping_grav0.1",
                         exeName = "P12_PBDPhysicsWithRotator.exe", otherArguments=rotatorArguments, recoverState=recoveryState )

    # generate 2 squishyBall colliding running parameters for AnkAPC01
    # recoveryState = r"Demo09_SquishyBall1_lowres\BowlBoundary\Hollow_2_MassSrping\run1\RecoveryStates\A00000900.json"
    # genRunningParameters("../../Bin/Run_AnkaPC01_2_squishyBall_MassSpring.cmd", "AnkaPC01", "Demo09_SquishyBall1_lowres\BowlBoundary\hollowMassSpring",
    #                      "Models_num_2_withCCD.json", "Parameters_num_2_withCCD.json", "run2_withCCD",
    #                      exeName = "P08_PBDPhysicsDemo_temp.exe",  recoverState=recoveryState )
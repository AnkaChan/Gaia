import json
import os
import subprocess
import time
from os.path import join
from pathlib import Path

machines = {
    "PC00":{
        "binaryFile": r"D:\Projects\Gaia_VBDCloth\Release\GAIA_VBDCloth.exe",
        "RepoPath": r"D:\Code\Graphics\Gaia",
        "OutputDirs": [
            r"D:\Data\VBD_cloth_Results",
        ]
    },

}

def genRunningParameters2(machineName, experimentPath, experimentName, modelsInfo, parametersInfo, modelsFile="Models.json", parametersFile="Parameters.json", outputDirId=0,
                         exeName = r"P10_EBDPhysicsDemo.exe", recoverState = None, relativeRecoverStatePath=True, otherArguments = None, writeCommand=True, runCommand=False, commandlineSeparator='\n'):
    machineInfo = machines[machineName]
    outputDir = machineInfo["OutputDirs"][outputDirId]

    parameterOutputDir = join(".", "Parameters", experimentPath, experimentName)
    os.makedirs(parameterOutputDir, exist_ok=True)

    json.dump(modelsInfo, open(join(parameterOutputDir, modelsFile), 'w'), indent=2)
    json.dump(parametersInfo, open(join(parameterOutputDir, parametersFile), 'w'), indent=2)

    modelsFileRemote = join(machineInfo["RepoPath"], "Simulator", "VBDCloth", "ParameterGen", "Parameters", experimentPath, experimentName, modelsFile)
    parametersFileRemote = join(machineInfo["RepoPath"], "Simulator", "VBDCloth", "ParameterGen", "Parameters", experimentPath, experimentName, parametersFile)

    cmd = [
        join(machineInfo["RepoPath"], "Binaries", exeName),
        modelsFileRemote,
        parametersFileRemote,
        join(outputDir, experimentPath, experimentName),
        '-R', machineInfo["RepoPath"]
    ]

    if recoverState is not None:
        if relativeRecoverStatePath:
            recoverState = join(outputDir, recoverState)
        cmd.extend(["-r", recoverState])

    if "temp" in exeName:
        print("Warning!!! temp exe is used!")


    if otherArguments is not None:
        if isinstance(otherArguments, list):
            cmd.extend(otherArguments)
        elif isinstance(otherArguments, str):
            cmd.append(otherArguments)
    elif exeName == r"P12_PBDPhysicsWithRotator.exe":
            print(exeName, " requires extra arguments!!!")

    if writeCommand:
        outputCmdFile = join(parameterOutputDir, "Run_" + machineName + ".cmd")
        f = open(outputCmdFile, 'w')
        f.write(" ".join(cmd))
        print("Command: ", commandlineSeparator.join(cmd))

        f.write("\ncmd /k")

    if runCommand:
        startT = time.time()
        subprocess.run(cmd)
        T = time.time() - startT
        print("Running command took: ", T)

    return cmd


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
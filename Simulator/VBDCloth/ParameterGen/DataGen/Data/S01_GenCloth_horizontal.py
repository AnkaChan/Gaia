import os

import numpy as np
import itertools
from S01_GenCloth import *

if __name__ == '__main__':
    # N = 5
    # size = (4, 4)
    # position = (0,0)
    #
    N = 223
    size = (100, 100)
    position = (0,0)

    # N =
    # N = 4
    # N=1000


    outFolder = "Data/SyntheticData"
    outObj = outFolder + f'/C{N}_horizontal.obj'
    outJson = outObj + ".vertexColoring.json"

    os.makedirs(outFolder, exist_ok=True)

    # size = (100, 150)

    # size = 10

    Y = np.linspace(-0.5 * size[0] + position[1], 0.5 * size[0] + position[1], N)
    X = np.linspace(-0.5 * size[1] + position[0], 0.5 * size[1] + position[0], N)

    U = np.linspace(0, 1, N, endpoint=True)
    V = np.linspace(0, 1, N, endpoint=True)

    X, Y = np.meshgrid(X, Y)
    U, V = np.meshgrid(U, V)

    Z = []
    for i in range(N):
        Z.append(np.zeros((N)))

    Z = np.array(Z)

    writeOBj_crossTesselation(outObj, N, X, Z, Y, U, V)

    writeOBj_crossTesselation(outObj, N, X, Z, Y, U, V)
    writeVertexColoring3(outJson, N)
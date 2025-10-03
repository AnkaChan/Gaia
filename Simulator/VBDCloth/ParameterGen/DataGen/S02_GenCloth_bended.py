import os

import numpy as np
import itertools
from S01_GenCloth import *

if __name__ == '__main__':
    N = 30
    # N = 4
    # N = 250

    outFolder = "Data/SyntheticData"
    outObj = outFolder + f'/C{N}_bended.obj'

    os.makedirs(outFolder, exist_ok=True)

    size = 100
    position = (0,0)

    X = []
    Y = []

    ts = np.linspace(0, 2*np.pi, N)
    for i in range(N):
        X.append(0.5 * size *np.sin(ts))
        Y.append(0.5 * size *np.cos(ts))

    X = np.array(X)
    Y = np.array(Y)
    # X = np.linspace(0.5 * size * s= + position[0], -0.5 * size + position[0], N)
    # Y = np.linspace(-0.5 * size + position[1], 0.5 * size + position[1], N)

    U = np.linspace(0, 1, N, endpoint=True)
    V = np.linspace(0, 1, N, endpoint=True)

    X, Y = np.meshgrid(X, Y)
    U, V = np.meshgrid(U, V)

    Z = []
    for i in range(N):
        Z.append(np.linspace(0, size, N))

    Z = np.array(Z)

    writeOBj_crossTesselation(outObj, N, X, Z, Y, U, V)
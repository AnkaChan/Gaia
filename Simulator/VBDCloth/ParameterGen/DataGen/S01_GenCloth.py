import os

import numpy as np
import itertools
import json
def writeOBj(outObj, N, X, Y, Z, U, V):
    file = open(outObj, 'w')
    for i, j in itertools.product(range(N), range(N)):
        file.write('v %f %f %f\n' % (X[i, j], Y[i, j], Z[i, j]))

    for i, j in itertools.product(range(N), range(N)):
        file.write('vt %f %f\n' % (U[i, j], V[i, j]))

    for i, j in itertools.product(range(0, N - 1), range(1, N)):
        vId = j + i * N
        file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + 1, vId + 1, vId + N + 1, vId + N + 1))
        file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + N + 1, vId + N + 1, vId + N, vId + N,))


def writeOBj_crossTesselation(outObj, N, X, Y, Z, U, V):
    file = open(outObj, 'w')
    for i, j in itertools.product(range(N), range(N)):
        file.write('v %f %f %f\n' % (X[i, j], Y[i, j], Z[i, j]))

    for i, j in itertools.product(range(N), range(N)):
        file.write('vt %f %f\n' % (U[i, j], V[i, j]))

    for i, j in itertools.product(range(0, N - 1), range(1, N)):
        vId = j + i * N

        if (j + i) % 2:
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + N + 1, vId + N + 1, vId + 1, vId + 1, ))
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + N, vId + N, vId + N + 1, vId + N + 1, ))
        else:
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + N, vId + N, vId + 1, vId + 1, ))
            file.write('f %d/%d %d/%d %d/%d\n' % (vId + N, vId + N, vId + N + 1, vId + N + 1, vId + 1, vId + 1, ))

def writeVertexColoring3(outJson, N):
    colors=[[],[],[]]
    index = 0
    for i in range(N):
        for j in range(N):
            if (i+j) % 2 == 0:
                if i % 2 == 0:
                    colors[0].append(index)
                else:
                    colors[1].append(index)
            else:
                colors[2].append(index)
            index += 1

    # save to outJson
    json.dump(colors, open(outJson, 'w'), indent=0)

if __name__ == '__main__':
    N = 50
    # N =
    N = 30
    # N=1000

    N = 100

    outFolder = "DataGen/Data/SyntheticData"
    outObj = outFolder + f'/C{N}.obj'
    outJson = outObj + ".vertexColoring.json"

    os.makedirs(outFolder, exist_ok=True)

    size = 100
    # size = 10
    position = (0,0)

    X = np.linspace(-0.5 * size + position[0], 0.5 * size + position[0], N)
    Y = np.linspace(-0.5 * size + position[1], 0.5 * size + position[1], N)

    U = np.linspace(0, 1, N, endpoint=True)
    V = np.linspace(0, 1, N, endpoint=True)

    X, Y = np.meshgrid(X, Y)
    U, V = np.meshgrid(U, V)

    Z = []
    for i in range(N):
        Z.append(np.linspace(0, size, N))

    Z = np.array(Z)

    writeOBj_crossTesselation(outObj, N, X, Z, Y, U, V)
    writeVertexColoring3(outJson, N)
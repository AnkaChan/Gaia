import subprocess
from pathlib import Path
from os.path import join
import json
import numpy as np
from copy import deepcopy
import itertools

def writeOBj_crossTesselation(outObj, res, X, Y, Z, U, V):
    file = open(outObj, 'w')
    for i, j in itertools.product(range(res[0]), range(res[1])):
        file.write('v %f %f %f\n' % (X[i, j], Y[i, j], Z[i, j]))

    for i, j in itertools.product(range(res[0]), range(res[1])):
        file.write('vt %f %f\n' % (U[i, j], V[i, j]))

    for i, j in itertools.product(range(0, res[0] - 1), range(1, res[1])):
        vId = j + i * res[1]

        if (j + i) % 2:
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + res[1] + 1, vId + res[1] + 1, vId + 1, vId + 1, ))
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + res[1], vId + res[1], vId + res[1] + 1, vId + res[1] + 1, ))
        else:
            file.write('f %d/%d %d/%d %d/%d\n' % (vId, vId, vId + res[1], vId + res[1], vId + 1, vId + 1, ))
            file.write('f %d/%d %d/%d %d/%d\n' % (vId + res[1], vId + res[1], vId + res[1] + 1, vId + res[1] + 1, vId + 1, vId + 1, ))

class ClothGenerator:
    def __init__(s):
        s.size = None        # [x,y]
        s.resolution = None  # [x,y]
        s.position = (0, 0)


    def gen(s):
        X = []
        Y = []
        for i in range(s.resolution[0]):
            # X.append(0.5 * s.size[0] * np.sin(i * 2 * np.pi / s.resolution[0]))
            X = np.linspace(-0.5 * s.size[0] + s.position[0], 0.5 * s.size[0] + s.position[0], s.resolution[0])

        for i in range(s.resolution[1]):
            # Y.append(0.5 * s.size[1] * np.cos(i * 2 * np.pi / s.resolution[1]))
            Y = np.linspace(-0.5 * s.size[1] + s.position[1], 0.5 * s.size[1] + s.position[1], s.resolution[1])

        X = np.array(X)
        Y = np.array(Y)
        # X = np.linspace(0.5 * size * s= + position[0], -0.5 * size + position[0], N)
        # Y = np.linspace(-0.5 * size + position[1], 0.5 * size + position[1], N)

        U = np.linspace(0, 1, s.resolution[0], endpoint=True)
        V = np.linspace(0, 1, s.resolution[1], endpoint=True)

        s.X, s.Y = np.meshgrid(Y, X)
        s.U, s.V = np.meshgrid(V, U)

        Z = []
        for i in range(s.resolution[0]):
            Z.append(np.linspace(0, 0, s.resolution[1]))

        s.Z = np.array(Z)

        # turn into faces and vertices lists
        s.verts = []
        s.uvs = []
        for i, j in itertools.product(range(s.resolution[0]), range(s.resolution[1])):
            s.verts.append(np.array([s.X[i, j], s.Y[i, j], s.Z[i, j]]))

        for i, j in itertools.product(range(s.resolution[0]), range(s.resolution[1])):
            s.uvs.append(np.array([s.U[i, j], s.V[i, j]]))

        s.faces = []
        res = s.resolution
        for i, j in itertools.product(range(0, res[0] - 1), range(1, res [1])):
            vId = j + i * res[1]

            if (j + i) % 2:
                s.faces.append(((vId, vId), (vId + res[1] + 1, vId + res[1] + 1), (vId + 1, vId + 1),))
                s.faces.append(((vId, vId), (vId + res[1], vId + res[1]), (vId + res[1] + 1, vId + res[1] + 1),))
            else:
                s.faces.append(((vId, vId), (vId + res[1], vId + res[1]), (vId + 1, vId + 1),))
                s.faces.append((((vId + res[1], vId + res[1]), (vId + res[1] + 1, vId + res[1] + 1), (vId + 1, vId + 1),)))

    def writeObj(s, outFile, ):
        with open(outFile, 'w+') as f:
            fp = Path(outFile)

            for i, v in enumerate(s.verts):

                if len(v) == 3:
                    f.write('v {:f} {:f} {:f}\n'.format(v[0], v[1], v[2]))
                elif len(v) == 6:
                    f.write('v {:f} {:f} {:f} {:f} {:f} {:f}\n'.format(v[0], v[1], v[2], v[3], v[4], v[5]))

            for vt in s.uvs:
                f.write('vt {:f} {:f}\n'.format(vt[0], vt[1]))


            for iF in range(len(s.faces)):
                # if facesToPreserve is not None and iF not in facesToPreserve:
                #     continue
                f.write('f')
                for fis in s.faces[iF]:
                    f.write(' {}'.format('/'.join([str(fi) for fi in fis])))
                f.write('\n')
            f.close()
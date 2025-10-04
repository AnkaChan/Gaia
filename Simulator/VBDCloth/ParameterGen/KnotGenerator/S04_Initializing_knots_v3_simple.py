import numpy as np
from scipy.spatial.transform import Rotation as R

from M01_Curves import *
from matplotlib import pyplot as plt
from M02_Geometry import *
# import polyscope as ps


def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = abs(x_limits[1] - x_limits[0])
    x_middle = np.mean(x_limits)
    y_range = abs(y_limits[1] - y_limits[0])
    y_middle = np.mean(y_limits)
    z_range = abs(z_limits[1] - z_limits[0])
    z_middle = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_middle - plot_radius, x_middle + plot_radius])
    ax.set_ylim3d([y_middle - plot_radius, y_middle + plot_radius])
    ax.set_zlim3d([z_middle - plot_radius, z_middle + plot_radius])
if __name__ == '__main__':

    knotControlPtsScale = 0.07
    # controlPtsScale = 1

    clothWidth = 0.045

    rodScale = [1, 1, 1]
    rodLengthScaleFactor = 1      # for increasing the resolution along the rod
    restposeLengthScaleFactor =1  # for making the total length shorter
    radiusScaleFactor = 1

    numberOfKnots = 3
    nCircles = 3
    circleStep = 0.12
    circleRadius = 0.09
    circleDivisions = 3

    resolution = (numberOfKnots* 500, 20)
    outFileName = f"knot_v3_{resolution[0]}x{resolution[1]}"

    headSize = 5
    # used to initialize the x-axis for the first coord sys
    seedVector = np.array([0, 0, 1 ])

    knotControlPts = []
    knotTranslations = [
        [3.8, 0, 0],
        [0.0, 0.5, 0 ],
        [-3.5, 0, 0],

    ]
    for i in range(numberOfKnots):
        pts = getOverhandKnotControlPts()
        shift = np.array(knotTranslations[i])
        pts = pts + shift
        knotControlPts.append(pts * knotControlPtsScale)

    controlPts = np.vstack(knotControlPts)

    # make the head and tail for pulling
    head = [
        [-0.5, -0.2, 0],
        [0, -0.1, 0],
        [0, -0.3, 0],
        [0.33, -0.2, 0.15]
    ]

    tail = [
        [-0.38, -0.25, -0.15],
        [0, -0.15, 0.1],
        [0.5, -0.15, 0.0],
    ]


    controlPts = np.vstack(
        [
            head,
            controlPts,
            tail
        ]
    )

    ts, tsSampled, lengths, coords, pts = getCurveLengthParameterization(controlPts, resolution[0], seedVector)

    print("X Min:", np.min(pts[:, 0]))
    print("X Max:", np.max(pts[:, 0]))

    ax = plt.axes(projection='3d')

    ax.plot3D(pts[:, 0], pts[:, 1], pts[:, 2], 'orange')
    colors = ['red', 'green', 'blue']
    # plot normals
    plotNormals = False
    if plotNormals:
        for i in range(1, resolution[0] - 1, 20):
            coord = coords[i - 1]
            for j in range(3):
                axis = coord[:, j]
                pt = pts[i, :]
                normalsVec = np.array([
                    pt,
                    pt + axis
                ])
                ax.plot3D(normalsVec[:, 0], normalsVec[:, 1], normalsVec[:, 2], colors[j])


    print("X Min:", np.min(pts[:, 0]))
    print("X Max:", np.max(pts[:, 0]))

    # visualize

    # ps.init()

    # generate some random nodes and edges between them
    edges = np.array([[i, i + 1] for i in range(len(pts) - 1)])

    # # visualize!
    # ps_net = ps.register_curve_network("my network", pts, edges)

    # # visualize some random vectors per-node
    # ps_net.add_vector_quantity("normal", coords[:,0], enabled=True)
    # ps_net.add_vector_quantity("tangent", coords[:,1], enabled=True)
    # ps_net.add_vector_quantity("binormal", coords[:,2], enabled=True)
    # ps.show()


    # ax = plt.axes(projection='3d')

    # ax.plot3D(pts[:, 0], pts[:, 1], pts[:, 2], 'orange')
    # colors = ['red', 'green', 'blue']
    # # plot normals
    # plotNormals = False
    # if plotNormals:
    #     for i in range(1, resolution - 1, 20):
    #         coord = coords[i - 1]
    #         for j in range(3):
    #             axis = coord[:, j]
    #             pt = pts[i, :]
    #             normalsVec = np.array([
    #                 pt,
    #                 pt + axis
    #             ])
    #             ax.plot3D(normalsVec[:, 0], normalsVec[:, 1], normalsVec[:, 2], colors[j])

    # set_axes_equal(ax)
    # ax.set_title('Curve')
    # plt.show()

    totalLength = lengths[-1]
    totalLengthScaled = totalLength * rodLengthScaleFactor

    clothGen = ClothGenerator()

    clothGen.size = (totalLength*1.01, clothWidth)
    clothGen.resolution = resolution
    clothGen.position = (clothGen.size[0]/2, clothGen.size[1]/2)

    clothGen.gen()
    vsStack = clothGen.verts
    # ps.init()
    # ps_mesh = ps.register_surface_mesh("my mesh", np.array(vsStack), np.array(clothGen.faces)[...,0]-1)
    # ps.show()
    vsStackNew = []
    stripeLength = clothGen.size[0]
    # deform vertices:
    startl = 0.01
    endl = 0.01
    origins = []
    new_ts = []
    ls = []
    v2s = []
    binormals = []
    for iV in range(len(vsStack)):
        v = np.array([vsStack[iV][2], vsStack[iV][1], vsStack[iV][0]])
        # v = vsStack[iV]
        l = v[1]
        l = startl + (totalLength - startl - endl) * l / stripeLength

        for iDim in range(3):
            if iDim == 1:
                # scale the length of the rod to be equal to the curve
                v[iDim] = l * restposeLengthScaleFactor
            else:
                v[iDim] = v[iDim] * rodScale[iDim] * radiusScaleFactor

        # convert length to parameters, and get coordinates system
        ls.append(l)
        t, coord = lenghtParameteriztioToCentripetalParameterization(l, tsSampled, lengths, coords)
        new_ts.append(t)
        origin, frame = computeCatmullRomCurveWithFrenetFrame(t, controlPts, ts)
        binormals.append(frame[:, 2])
        origins.append(origin)
        vDeformed = origin + v[0] * frame[:, 0] + v[2] * frame[:, 2]
        v2s.append(v[2])
        vsStackNew.append(vDeformed.tolist())

    # write the original vertex positions
    origins = np.array(origins)
    center = origins[::resolution[1]]
    new_ts = np.array(new_ts)[::resolution[1]]
    ls = np.array(ls)[::resolution[1]]

    diff = center[1:] - center[:-1]
    diff_norm = np.linalg.norm(diff, axis=1)
    diff = diff / diff_norm[:, None]

    tangent = np.zeros((len(center), 3))
    tangent[1:-1] = diff[:-1] + diff[1:]
    tangent[0] = diff[0]
    tangent[-1] = diff[-1]
    tangent_norm = np.linalg.norm(tangent, axis=1)
    tangent = tangent / tangent_norm[:, None]

    normal = np.zeros((len(center), 3))
    normal[1:-1] = diff[1:] - diff[:-1]
    normal[0] = normal[1]
    normal[-1] = normal[-2]
    normal_norm = np.linalg.norm(normal, axis=1)
    eps = 1e-6
    mask = normal_norm < eps
    flat =(np.where(mask))[0]
    normal[~mask] = normal[~mask] / normal_norm[~mask, None]
    for i in flat:
        normal[i] = normal[i-1]

    for i in range(1, len(normal)):
        if np.dot(normal[i], normal[i-1]) < 0:
            normal[i] = -normal[i]

    binormals = np.cross(tangent, normal)
    # ps.init()
    # ps_net = ps.register_curve_network("my network", np.array(origins), np.array([[i, i + 1] for i in range(len(origins) - 1)]))
    # ps_net.add_vector_quantity("normal", normal, enabled=True)
    # ps_net.add_vector_quantity("binormal", binormals, enabled=True)
    # ps.show()
    v2s = np.array(v2s)
    binormals = np.repeat(binormals, resolution[1], axis=0)

    edgeLensOrg = []
    for f in clothGen.faces:
        for iFV in range(len(f)):
            fv1 = f[iFV][0]-1
            fv2 = f[(iFV+1) %3 ][0]-1
            l = np.linalg.norm(clothGen.verts[fv1] - clothGen.verts[fv2])
            edgeLensOrg.append(l)

    vsStackNew = origins + binormals * v2s[:, None]

    clothGen.writeObj(join("Data", outFileName + "_restShape.obj"))
    clothGen.verts = vsStackNew
    edgeLensDeformed = []
    for f in clothGen.faces:
        for iFV in range(len(f)):
            fv1 = f[iFV][0]-1
            fv2 = f[(iFV+1) %3 ][0]-1
            l = np.linalg.norm(clothGen.verts[fv1] - clothGen.verts[fv2])
            edgeLensDeformed.append(l)

    print((np.array(edgeLensDeformed)/np.array(edgeLensOrg)).mean())

    clothGen.writeObj(join("Data", outFileName + ".obj"))


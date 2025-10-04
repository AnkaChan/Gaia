import numpy as np
from scipy.spatial.transform import Rotation as R

from M01_Curves import *
from matplotlib import pyplot as plt
from M02_Geometry import *



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



    outFile = r'..\Data\mesh_models\ModelsDifferentRes\knots\OverhandKnot_fractal3.t'
    outColoringFile = outFile  + '.coloring.json'

    controlPtsScale = 0.025
    # controlPtsScale = 1

    rodScale = [1.8, 1, 1.8]
    rodLengthScaleFactor = 1.3      # for increasing the resolution along the rod
    restposeLengthScaleFactor = 0.7 # for making the total length shorter
    radiusScaleFactor = 1.35

    numberOfKnots = 6
    nCircles = 3
    circleStep = 0.1
    circleRadius = 0.05
    circleDivisions = 8

    resolution = numberOfKnots * 300

    headSize = 5
    shift = np.array((-3.5, 0, 0))
    # used to initialize the x-axis for the first coord sys
    seedVector = np.array([0, 0, 1 ])

    controlPts = []
    for i in range(numberOfKnots):
        pts = getOverhandKnotControlPts()
        pts = pts + shift * i
        controlPts.append(pts * controlPtsScale)

    controlPts = np.vstack(controlPts)

    # make the head and tail for pulling
    head = [[-1.2, -0.3, 0],
            [-0.9, -0.3, 0]
            ]

    for iCirles in range(1, nCircles+1):
        for iDiv in range(0, circleDivisions):
            head.append(
                [head[1][0] + (iCirles + iDiv / circleDivisions) * circleStep,
                 head[1][1] + circleRadius * np.cos(2 * np.pi * iDiv / circleDivisions),
                 head[1][2] + circleRadius * np.sin(2 * np.pi * iDiv / circleDivisions)]
            )
        # head.append(
        #     [head[1][0] + iCirles * circleStep, head[1][1] + circleRadius,   head[1][2] ]
        # )
        # head.append(
        #     [head[1][0] + (iCirles + 0.25) * circleStep, head[1][1],   head[1][2] + circleRadius ]
        # )
        # head.append(
        #     [head[1][0] + (iCirles + 0.5) * circleStep, head[1][1] - circleRadius,   head[1][2] ]
        # )
        # head.append(
        #     [head[1][0] + (iCirles + 0.75) * circleStep, head[1][1], head[1][2] - circleRadius]
        # )

    head.append( [-0.3, -0.3, 0.0])
    head.append( [0.1, -0.3, 0.0])
    head.append( [0.1, -0.3, 0.1])
    head.append( [0.1, -0.1, 0.1])

    head = np.array(head)
    # head[:, 0] = head[:, 0] + headSize * controlPtsScale
    tail = [
        [-0.8, -0.3, 0.1],
        [-0.8, -0.3, 0.0],
        [-0.7, -0.3, 0.0],
        [-0.4, -0.3, 0.0],
    ]

    end = tail[-1]

    for iCirles in range(1, nCircles+1):
        tail.append(
            [end[0] + iCirles * circleStep, end[1] + circleRadius,   end[2] ]
        )
        tail.append(
            [end[0] + (iCirles + 0.25) * circleStep, end[1],   end[2] + circleRadius ]
        )
        tail.append(
            [end[0] + (iCirles + 0.5) * circleStep, end[1] - circleRadius,   end[2] ]
        )
        tail.append(
            [end[0] + (iCirles + 0.75) * circleStep, end[1], end[2] - circleRadius]
        )

    tail.append( [0.2, -0.3, 0.0])
    tail.append( [0.3, -0.3, 0.0])
    # tail[:, 0] = tail[:, 0] - headSize * controlPtsScale
    tail = np.array(tail)

    controlPts = np.vstack(
        [
            head,
            controlPts,
            tail
        ]
    )

    ts, tsSampled, lengths, coords, pts = getCurveLengthParameterization(controlPts, resolution, seedVector)

    print("X Min:", np.min(pts[:, 0]))
    print("X Max:", np.max(pts[:, 0]))

    ax = plt.axes(projection='3d')

    ax.plot3D(pts[:, 0], pts[:, 1], pts[:, 2], 'orange')
    colors = ['red', 'green', 'blue']
    # plot normals
    plotNormals = False
    if plotNormals:
        for i in range(1, resolution - 1, 20):
            coord = coords[i - 1]
            for j in range(3):
                axis = coord[:, j]
                pt = pts[i, :]
                normalsVec = np.array([
                    pt,
                    pt + axis
                ])
                ax.plot3D(normalsVec[:, 0], normalsVec[:, 1], normalsVec[:, 2], colors[j])

    set_axes_equal(ax)
    ax.set_title('Curve')
    plt.show()

    totalLength = lengths[-1]
    totalLengthScaled = totalLength * rodLengthScaleFactor

    numOfStack = int(np.ceil(totalLengthScaled)) - 1

    clothGen = ClothGenerator()

    clothGen.size = (8, 0.005)
    clothGen.resolution = (1000, 10)
    clothGen.position = (clothGen.size[0]/2, clothGen.size[1]/2)

    clothGen.gen()
    vsStack = clothGen.verts

    vsStackNew = []
    stripeLength = clothGen.size[0]
    # deform vertices:
    for iV in range(len(vsStack)):
        # v = np.array([vsStack[iV][0], vsStack[iV][2], -vsStack[iV][1]])
        v = vsStack[iV]
        l = v[1]
        l = totalLength * l / stripeLength

        for iDim in range(3):
            if iDim == 1:
                # scale the length of the rod to be equal to the curve
                v[iDim] = l * restposeLengthScaleFactor
            else:
                v[iDim] = v[iDim] * rodScale[iDim] * radiusScaleFactor

        # convert length to parameters, and get coordinates system
        t, coord = lenghtParameteriztioToCentripetalParameterization(l, tsSampled, lengths, coords)

        origin = computeCatmullRomCurve(t, controlPts, ts)

        vDeformed = origin + v[0] * coord[:, 0] + v[2] * coord[:, 2]
        vsStackNew.append(vDeformed.tolist())

    # write the original vertex positions

    clothGen.verts = vsStackNew

    clothGen.writeObj("Data/knot.obj")



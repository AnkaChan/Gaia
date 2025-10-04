import numpy as np
from scipy.spatial.transform import Rotation as R

def interpolate(x1, x2, t, t0, t1):
    if t0 == t1:
        return x1

    return x1 * (t1-t)/(t1-t0) + x2 * (t-t0)/(t1-t0)

def cubicCatmullRom(t, ps, ts):
    L01 = interpolate(ps[0], ps[1], t, ts[0], ts[1])
    L12 = interpolate(ps[1], ps[2], t, ts[1], ts[2])
    L23 = interpolate(ps[2], ps[3], t, ts[2], ts[3])

    L012 = interpolate(L01, L12, t, ts[0], ts[2])
    L123 = interpolate(L12, L23, t, ts[1], ts[3])

    C12 = interpolate(L012, L123, t, ts[1], ts[2])

    return C12

def computeCatmullRomCurve(t, controlPoints, ts):
    a = np.max(np.where(t >= ts))
    b = np.min(np.where(t <= ts))

    if a == b:
        return controlPoints[a]
    assert a == b-1

    if a - 1 < 0:
        selectedControlPts = controlPoints[0:b + 2]
        selectedControlPts = np.vstack([selectedControlPts[0:1], selectedControlPts])
        tsSelected = ts[0:b + 2].tolist()
        tsSelected.insert(0, tsSelected[0])
    elif b +1 >= controlPoints.shape[0]:
        selectedControlPts = controlPoints[a - 1:]
        selectedControlPts = np.vstack([selectedControlPts, selectedControlPts[-1:]])
        tsSelected = ts[a-1:].tolist()
        tsSelected.insert(-1, tsSelected[-1])
    else:
        selectedControlPts = controlPoints[a - 1:b + 2]
        tsSelected = ts[a-1:b+2]

    return cubicCatmullRom(t, selectedControlPts, tsSelected)

def computeCatmullRomCurveWithFrenetFrame(t, controlPoints, ts):
    a = np.max(np.where(t >= ts))
    b = np.min(np.where(t <= ts))

    if a == b:
        b = a + 1
    assert a == b-1

    if a - 1 < 0:
        selectedControlPts = controlPoints[0:b + 2]
        selectedControlPts = np.vstack([selectedControlPts[0:1], selectedControlPts])
        tsSelected = ts[0:b + 2].tolist()
        tsSelected.insert(0, tsSelected[0])
    elif b +1 >= controlPoints.shape[0]:
        selectedControlPts = controlPoints[a - 1:]
        selectedControlPts = np.vstack([selectedControlPts, selectedControlPts[-1:]])
        tsSelected = ts[a-1:].tolist()
        tsSelected.insert(-1, tsSelected[-1])
    else:
        selectedControlPts = controlPoints[a - 1:b + 2]
        tsSelected = ts[a-1:b+2]

    x = cubicCatmullRom(t, selectedControlPts, tsSelected)

    t01 = tsSelected[0]
    t02 = tsSelected[1]
    t03 = tsSelected[2]
    t04 = tsSelected[3]

    p11 = selectedControlPts[0][0]
    p12 = selectedControlPts[0][1]
    p13 = selectedControlPts[0][2]

    p21 = selectedControlPts[1][0]
    p22 = selectedControlPts[1][1]
    p23 = selectedControlPts[1][2]

    p31 = selectedControlPts[2][0]
    p32 = selectedControlPts[2][1]
    p33 = selectedControlPts[2][2]

    p41 = selectedControlPts[3][0]
    p42 = selectedControlPts[3][1]
    p43 = selectedControlPts[3][2]

    tangent_n = [
        (((((p11 * (t - t02)) / (t01 - t02) - (p21 * (t - t01)) / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                ((p21 * (t - t03)) / (t02 - t03) - (p31 * (t - t02)) / (t02 - t03)) * (t - t01)) / (t01 - t03)) / (
                 t02 - t03) - (
                 (((p21 * (t - t03)) / (t02 - t03) - (p31 * (t - t02)) / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                     ((p31 * (t - t04)) / (t03 - t04) - (p41 * (t - t03)) / (t03 - t04)) * (t - t02)) / (t02 - t04)) / (
                 t02 - t03) + ((t - t03) * (
                ((p11 * (t - t02)) / (t01 - t02) - (p21 * (t - t01)) / (t01 - t02)) / (t01 - t03) - (
                    (p21 * (t - t03)) / (t02 - t03) - (p31 * (t - t02)) / (t02 - t03)) / (t01 - t03) + (
                            (p11 / (t01 - t02) - p21 / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                            (p21 / (t02 - t03) - p31 / (t02 - t03)) * (t - t01)) / (t01 - t03))) / (t02 - t03) - (
                 (t - t02) * (((p21 * (t - t03)) / (t02 - t03) - (p31 * (t - t02)) / (t02 - t03)) / (t02 - t04) - (
                     (p31 * (t - t04)) / (t03 - t04) - (p41 * (t - t03)) / (t03 - t04)) / (t02 - t04) + (
                                          (p21 / (t02 - t03) - p31 / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                                          (p31 / (t03 - t04) - p41 / (t03 - t04)) * (t - t02)) / (t02 - t04))) / (
                 t02 - t03)),
        (((((p12 * (t - t02)) / (t01 - t02) - (p22 * (t - t01)) / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                ((p22 * (t - t03)) / (t02 - t03) - (p32 * (t - t02)) / (t02 - t03)) * (t - t01)) / (t01 - t03)) / (
                 t02 - t03) - (
                 (((p22 * (t - t03)) / (t02 - t03) - (p32 * (t - t02)) / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                     ((p32 * (t - t04)) / (t03 - t04) - (p42 * (t - t03)) / (t03 - t04)) * (t - t02)) / (t02 - t04)) / (
                 t02 - t03) + ((t - t03) * (
                ((p12 * (t - t02)) / (t01 - t02) - (p22 * (t - t01)) / (t01 - t02)) / (t01 - t03) - (
                    (p22 * (t - t03)) / (t02 - t03) - (p32 * (t - t02)) / (t02 - t03)) / (t01 - t03) + (
                            (p12 / (t01 - t02) - p22 / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                            (p22 / (t02 - t03) - p32 / (t02 - t03)) * (t - t01)) / (t01 - t03))) / (t02 - t03) - (
                 (t - t02) * (((p22 * (t - t03)) / (t02 - t03) - (p32 * (t - t02)) / (t02 - t03)) / (t02 - t04) - (
                     (p32 * (t - t04)) / (t03 - t04) - (p42 * (t - t03)) / (t03 - t04)) / (t02 - t04) + (
                                          (p22 / (t02 - t03) - p32 / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                                          (p32 / (t03 - t04) - p42 / (t03 - t04)) * (t - t02)) / (t02 - t04))) / (
                 t02 - t03)),
        (((((p13 * (t - t02)) / (t01 - t02) - (p23 * (t - t01)) / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                ((p23 * (t - t03)) / (t02 - t03) - (p33 * (t - t02)) / (t02 - t03)) * (t - t01)) / (t01 - t03)) / (
                 t02 - t03) - (
                 (((p23 * (t - t03)) / (t02 - t03) - (p33 * (t - t02)) / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                     ((p33 * (t - t04)) / (t03 - t04) - (p43 * (t - t03)) / (t03 - t04)) * (t - t02)) / (t02 - t04)) / (
                 t02 - t03) + ((t - t03) * (
                ((p13 * (t - t02)) / (t01 - t02) - (p23 * (t - t01)) / (t01 - t02)) / (t01 - t03) - (
                    (p23 * (t - t03)) / (t02 - t03) - (p33 * (t - t02)) / (t02 - t03)) / (t01 - t03) + (
                            (p13 / (t01 - t02) - p23 / (t01 - t02)) * (t - t03)) / (t01 - t03) - (
                            (p23 / (t02 - t03) - p33 / (t02 - t03)) * (t - t01)) / (t01 - t03))) / (t02 - t03) - (
                 (t - t02) * (((p23 * (t - t03)) / (t02 - t03) - (p33 * (t - t02)) / (t02 - t03)) / (t02 - t04) - (
                     (p33 * (t - t04)) / (t03 - t04) - (p43 * (t - t03)) / (t03 - t04)) / (t02 - t04) + (
                                          (p23 / (t02 - t03) - p33 / (t02 - t03)) * (t - t04)) / (t02 - t04) - (
                                          (p33 / (t03 - t04) - p43 / (t03 - t04)) * (t - t02)) / (t02 - t04))) / (
                 t02 - t03))
        ]

    tangent_n = np.array(tangent_n)
    tangent_norm =  np.linalg.norm(tangent_n)
    tangent_n = tangent_n / tangent_norm

    normal = [
        -((
                      2 * p11 * t02 ** 2 - 2 * p11 * t03 ** 2 - 2 * p21 * t01 ** 2 + 2 * p21 * t03 ** 2 + 2 * p31 * t01 ** 2 - 2 * p31 * t02 ** 2 - 4 * p11 * t * t02 + 4 * p11 * t * t03 + 4 * p21 * t * t01 - 4 * p21 * t * t03 - 4 * p31 * t * t01 + 4 * p31 * t * t02) / (
                      (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) - (
                      2 * p21 * t03 ** 2 - 2 * p21 * t04 ** 2 - 2 * p31 * t02 ** 2 + 2 * p31 * t04 ** 2 + 2 * p41 * t02 ** 2 - 2 * p41 * t03 ** 2 - 4 * p21 * t * t03 + 4 * p21 * t * t04 + 4 * p31 * t * t02 - 4 * p31 * t * t04 - 4 * p41 * t * t02 + 4 * p41 * t * t03) / (
                      (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04)) - (
                      2 * (t - t03) * (p11 * t02 - p11 * t03 - p21 * t01 + p21 * t03 + p31 * t01 - p31 * t02)) / (
                      (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) + (
                      2 * (t - t02) * (p21 * t03 - p21 * t04 - p31 * t02 + p31 * t04 + p41 * t02 - p41 * t03)) / (
                      (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04))) / tangent_norm,
        - ((
                       2 * p12 * t02 ** 2 - 2 * p12 * t03 ** 2 - 2 * p22 * t01 ** 2 + 2 * p22 * t03 ** 2 + 2 * p32 * t01 ** 2 - 2 * p32 * t02 ** 2 - 4 * p12 * t * t02 + 4 * p12 * t * t03 + 4 * p22 * t * t01 - 4 * p22 * t * t03 - 4 * p32 * t * t01 + 4 * p32 * t * t02) / (
                       (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) - (
                       2 * p22 * t03 ** 2 - 2 * p22 * t04 ** 2 - 2 * p32 * t02 ** 2 + 2 * p32 * t04 ** 2 + 2 * p42 * t02 ** 2 - 2 * p42 * t03 ** 2 - 4 * p22 * t * t03 + 4 * p22 * t * t04 + 4 * p32 * t * t02 - 4 * p32 * t * t04 - 4 * p42 * t * t02 + 4 * p42 * t * t03) / (
                       (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04)) - (
                       2 * (t - t03) * (p12 * t02 - p12 * t03 - p22 * t01 + p22 * t03 + p32 * t01 - p32 * t02)) / (
                       (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) + (
                       2 * (t - t02) * (p22 * t03 - p22 * t04 - p32 * t02 + p32 * t04 + p42 * t02 - p42 * t03)) / (
                       (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04))) / tangent_norm,
        - ((
                       2 * p13 * t02 ** 2 - 2 * p13 * t03 ** 2 - 2 * p23 * t01 ** 2 + 2 * p23 * t03 ** 2 + 2 * p33 * t01 ** 2 - 2 * p33 * t02 ** 2 - 4 * p13 * t * t02 + 4 * p13 * t * t03 + 4 * p23 * t * t01 - 4 * p23 * t * t03 - 4 * p33 * t * t01 + 4 * p33 * t * t02) / (
                       (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) - (
                       2 * p23 * t03 ** 2 - 2 * p23 * t04 ** 2 - 2 * p33 * t02 ** 2 + 2 * p33 * t04 ** 2 + 2 * p43 * t02 ** 2 - 2 * p43 * t03 ** 2 - 4 * p23 * t * t03 + 4 * p23 * t * t04 + 4 * p33 * t * t02 - 4 * p33 * t * t04 - 4 * p43 * t * t02 + 4 * p43 * t * t03) / (
                       (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04)) - (
                       2 * (t - t03) * (p13 * t02 - p13 * t03 - p23 * t01 + p23 * t03 + p33 * t01 - p33 * t02)) / (
                       (t01 - t02) * (t01 - t03) * (t02 - t03) ** 2) + (
                       2 * (t - t02) * (p23 * t03 - p23 * t04 - p33 * t02 + p33 * t04 + p43 * t02 - p43 * t03)) / (
                       (t02 - t03) ** 2 * (t02 - t04) * (t03 - t04))) / tangent_norm,
    ]

    normal = np.array(normal)

    normal =  normal - tangent_n * np.dot(tangent_n, normal) / (tangent_norm*tangent_norm)

    normal = normal / np.linalg.norm(normal)

    binoraml = np.cross(tangent_n, normal)

    binoraml = binoraml / np.linalg.norm(binoraml)

    return x, np.hstack([tangent_n[...,None], normal[...,None], binoraml[..., None]])
def centripetalparameteriztion(controlPoints):
    ts = [0]

    for i in range(1, controlPoints.shape[0]):
        ts.append(ts[i-1] + np.sqrt(np.linalg.norm(controlPoints[i-1] - controlPoints[i])))

    return np.array(ts)

def computeDirection(t, controlPoints, ts, dt):
    if t - dt >=0:
        diff = computeCatmullRomCurve(t + dt, controlPoints, ts) - computeCatmullRomCurve(t - dt, controlPoints, ts)
    else:
        diff = computeCatmullRomCurve(t + dt, controlPoints, ts) - computeCatmullRomCurve(t, controlPoints, ts)

    return diff / np.linalg.norm(diff)


def getCurveLengthParameterization(controlPts, resolution, seedVector, dtResolution= 10):
    ts = centripetalparameteriztion(controlPts)

    tMax = ts[-1]

    pts = []

    lengths = [0]

    totalLength = 0
    for i in range(resolution):
        t = tMax * (i / resolution)

        p = computeCatmullRomCurve(t, controlPts, ts)
        pts.append(p)

        if i > 0:
            length = np.linalg.norm(pts[-1] - pts[-2])
            totalLength += length
            lengths.append(totalLength)
    print("Max t: ", tMax)
    print("Total length: ", totalLength)

    # compute normals and coordinates system at each
    normals = []
    dt = tMax / (resolution * dtResolution)

    coords = []
    tsSampled = []
    for i in range(resolution):
        t = tMax * (i / resolution)
        tsSampled.append(t)
        if t + dt < tMax:
            normal = computeDirection(t, controlPts, ts, dt)
            normals.append(normal)

            if len(coords) == 0:
                yAxis = normal
                # orthonormalize x-axis
                xAxis = seedVector - yAxis * (np.dot(seedVector, yAxis))
                xAxis = xAxis / np.linalg.norm(xAxis)
                zAxis = np.cross(xAxis, yAxis)
                newCoord = np.array([
                    xAxis,
                    yAxis,
                    zAxis
                ]).transpose()

            else:
                coord = coords[-1]
                newYAxis = normal
                oldYAxis = coord[:, 1]
                # compute new coord by rotating
                angle = np.arccos(np.dot(oldYAxis, newYAxis))
                if angle > 0.01:
                    axis = np.cross(oldYAxis, newYAxis)
                    axis = axis / np.linalg.norm(axis)

                    r = R.from_rotvec(angle * axis)
                    rMat = r.as_matrix()

                    newCoord = rMat @ coord
                    for iAxis in range(3):
                        newCoord[:, iAxis] = newCoord[:, iAxis] / np.linalg.norm(newCoord[:, iAxis])
                else:
                    newCoord = coord

            assert np.linalg.norm(normal - newCoord[:, 1]) < 1e-2

            r = R.from_matrix(newCoord)
            coords.append(r.as_matrix())

    return np.array(ts), np.array(tsSampled), np.array(lengths), np.array(coords), np.array(pts)


def getCurveLengthParameterizationWithFrenetFrame(controlPts, resolution, seedVector, dtResolution= 10):
    ts = centripetalparameteriztion(controlPts)

    tMax = ts[-1]

    pts = []

    lengths = [0]

    totalLength = 0
    coords = []

    for i in range(resolution):
        t = tMax * (i / resolution)

        p, coord = computeCatmullRomCurveWithFrenetFrame(t, controlPts, ts)
        coords.append(coord)
        pts.append(p)

        if i > 0:
            length = np.linalg.norm(pts[-1] - pts[-2])
            totalLength += length
            lengths.append(totalLength)
    print("Max t: ", tMax)
    print("Total length: ", totalLength)

    # compute normals and coordinates system at each
    normals = []
    dt = tMax / (resolution * dtResolution)

    tsSampled = []
    for i in range(resolution):
        t = tMax * (i / resolution)
        tsSampled.append(t)
        # if t + dt < tMax:
        #     normal = computeDirection(t, controlPts, ts, dt)
        #     normals.append(normal)
        #
        #     if len(coords) == 0:
        #         yAxis = normal
        #         # orthonormalize x-axis
        #         xAxis = seedVector - yAxis * (np.dot(seedVector, yAxis))
        #         xAxis = xAxis / np.linalg.norm(xAxis)
        #         zAxis = np.cross(xAxis, yAxis)
        #         newCoord = np.array([
        #             xAxis,
        #             yAxis,
        #             zAxis
        #         ]).transpose()
        #
        #     else:
        #         coord = coords[-1]
        #         newYAxis = normal
        #         oldYAxis = coord[:, 1]
        #         # compute new coord by rotating
        #         angle = np.arccos(np.dot(oldYAxis, newYAxis))
        #         if angle > 0.01:
        #             axis = np.cross(oldYAxis, newYAxis)
        #             axis = axis / np.linalg.norm(axis)
        #
        #             r = R.from_rotvec(angle * axis)
        #             rMat = r.as_matrix()
        #
        #             newCoord = rMat @ coord
        #             for iAxis in range(3):
        #                 newCoord[:, iAxis] = newCoord[:, iAxis] / np.linalg.norm(newCoord[:, iAxis])
        #         else:
        #             newCoord = coord
        #
        #     assert np.linalg.norm(normal - newCoord[:, 1]) < 1e-2
        #
        #     r = R.from_matrix(newCoord)
        #     coords.append(r.as_matrix())

    return np.array(ts), np.array(tsSampled), np.array(lengths), np.array(coords), np.array(pts)


def lenghtParameteriztioToCentripetalParameterization(l, tsSampled, lengths, coords):
    a = np.max(np.where(l >= lengths))
    b = np.min(np.where(l <= lengths))

    if a == b:
        t = tsSampled[a]
    else:
        assert a == b-1
        t = tsSampled[a] + (tsSampled[b] - tsSampled[a]) * (l - lengths[a]) / (lengths[b] - lengths[a])

    return t, coords[a]

def getOverhandKnotControlPts():
    pts = [
        # [1.4, 0, 0],
        [1, 0, 0],
        [0, 1, 0.0],
        [-1, 0, 0.0],
        [0, -1, 0.0],
        [1, 0, 1.5],
        [1, 1, 1],
        [1, 1, -0.75],
        [0, 0, -0.75],
        [0, 0, 1.25],
        [-1, 0, 1.5],
        [-2, 0, 0],
        # [-3, 0, 0],
        # [-4, 0, 0],

    ]

    return pts

if __name__ == '__main__':
    resolution = 500

    controlPts = np.array([
        [1, 0],
        [1, 2],
        [2, 2],
        [2, 1],
    ])

    ts = centripetalparameteriztion(controlPts)

    print("Parameterization: ", ts)

    tMax = ts[-1]

    for i in range(resolution):
        resolution


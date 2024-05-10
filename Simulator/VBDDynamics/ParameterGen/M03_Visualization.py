import pyvista as pv
import numpy as np
import os
from pyvista import examples
from tqdm import tqdm


def getBoxMesh(box):
    grid = examples.cells.Voxel()
    points = []
    for i in range(8):
        points.append([box[i % 2][0], box[i // 4][1], box[(i // 2) % 2][2]])
    points = np.array(points)
    grid.points = points
    return grid


def drawStatic(shape, window_size, func, args=[], save_path=None, show=True, cameraParams=None):
    n_row = shape[0]
    n_col = shape[1]
    n_size = n_row * n_col
    p = pv.Plotter(shape=shape, window_size=window_size)

    cameras = []
    for i, arg in enumerate(args):
        if i >= n_size:
            break
        row = i // n_col
        col = i % n_col
        p.subplot(row, col)
        if cameraParams is not None:
            camera = cameraParams[i]
            if camera is not None:
                p.camera.position = camera["position"]
                p.camera.focal_point = camera["focal_point"]
                p.camera.up = camera["up"]

        func(p, i, arg)
        cameras.append(p.camera)
    if save_path is not None:
        dirpath = os.path.dirname(os.path.abspath(save_path))
        if not os.path.exists(dirpath):
            os.makedirs(dirpath)
        p.export_html(save_path)
    if show:
        p.show()
    return cameras


def drawMovie(
    shape,
    window_size,
    save_path,
    dynamic_func,
    end_frame,
    start_frame=0,
    args=[],
    cameras=[],
    framerate=60,
    quality=5,
    static_func=None,
    showFrame=True,
    showPreview=False,
):
    dirpath = os.path.dirname(os.path.abspath(save_path))
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)
    n_row = shape[0]
    n_col = shape[1]
    n_size = n_row * n_col
    p = pv.Plotter(shape=shape, window_size=window_size, off_screen=not showPreview)
    p.open_movie(
        save_path,
        framerate=framerate,
        quality=quality,
    )
    for i, arg in enumerate(args):
        if i >= n_size:
            break
        row = i // n_col
        col = i % n_col
        p.subplot(row, col)
        static_func(p, i, arg)
    for frame in tqdm(range(start_frame, end_frame), desc="Creating MP4..."):
        actors = []
        for i, arg in enumerate(args):
            if i >= n_size:
                break
            row = i // n_col
            col = i % n_col
            p.subplot(row, col)
            new_actors = dynamic_func(p, frame, i, arg)
            if new_actors is not None:
                actors += new_actors
            if i < len(cameras):
                p.camera = cameras[i]
        if showFrame:
            frame_text = p.add_text(
                f"{frame:04d}", position="lower_right", font_size=10
            )
        p.write_frame()
        p.remove_actor(actors)
        if showFrame:
            p.remove_actor(frame_text)
    p.close()


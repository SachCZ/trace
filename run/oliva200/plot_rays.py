import rayvis
from matplotlib import pyplot as plt
import sys
import os
import msgpack
import numpy as np


def plot_rays(folder):
    fig, axes = plt.subplots()
    with open(os.path.join(folder, "output/mesh.mfem")) as f:
        mesh = rayvis.read_mfem_mesh(f)
    with open(os.path.join(folder, "output/power.gf")) as f:
        grid_function = rayvis.read_grid_function(f, mesh)
    with open(os.path.join(folder, "output/rays.msgpack"), "rb") as f:
        rays = rayvis.read_msgpack_rays(f)

    poly_collection = rayvis.plot_grid_function(axes, grid_function, cmap="GnBu")
    fig.colorbar(poly_collection)
    rayvis.plot_mesh(axes, mesh, linewidth=0.5)
    rayvis.plot_rays(axes, rays, linewidth=0.5, color="red")

    plt.show()
    plt.savefig(os.path.join(folder, "output/rays.png"))


if __name__ == '__main__':
    plot_rays(".")

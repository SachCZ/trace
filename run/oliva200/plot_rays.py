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
        mesh.nodes = mesh.nodes * 1e4
    with open(os.path.join(folder, "output/power.gf")) as f:
        grid_function = rayvis.read_grid_function(f, mesh)
        grid_function.values = np.abs(grid_function.values)
    with open(os.path.join(folder, "output/rays.msgpack"), "rb") as f:
        rays = rayvis.read_msgpack_rays(f)
        for ray in rays:
            ray.x = np.asarray(ray.x) * 1e4
            ray.y = np.asarray(ray.y) * 1e4

    poly_collection = rayvis.plot_grid_function(axes, grid_function)
    fig.colorbar(poly_collection)
    rayvis.plot_mesh(axes, mesh, linewidth=0.5)
    rayvis.plot_rays(axes, rays[::40], linewidth=0.8, color="red")
    plt.xlabel("x[μm]")
    plt.ylabel("y[μm]")
    plt.title("P [J$\\cdot$s$^{-1}$]")

    plt.savefig(os.path.join(folder, "output/rays.png"))
    plt.show()


if __name__ == '__main__':
    plot_rays(".")

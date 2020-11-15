import rayvis
from matplotlib import pyplot as plt
import sys
import os


def get_rays_subset(rays, n):
    each_n = round(len(rays) / n)
    return (ray for i, ray in enumerate(rays) if i % each_n == 0)


def plot_rays(folder, rays_count):
    fig, axes = plt.subplots()
    with open(os.path.join(folder, "input/mesh.mesh")) as f:
        mesh = rayvis.read_mfem_mesh(f)
    with open(os.path.join(folder, "output/absorbed_energy.gf")) as f:
        grid_function = rayvis.read_grid_function(f, mesh)
    with open(os.path.join(folder, "output/rays.msgpack"), "rb") as f:
        rays = rayvis.read_msgpack_rays(f)

    poly_collection = rayvis.plot_grid_function(axes, grid_function, cmap="GnBu")
    fig.colorbar(poly_collection)
    rayvis.plot_mesh(axes, mesh, linewidth=0.5)
    rayvis.plot_rays(axes, get_rays_subset(rays, rays_count), linewidth=0.5, color="red")

    plt.show()


if __name__ == '__main__':
    plot_rays(sys.argv[1], int(sys.argv[2]))

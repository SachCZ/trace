import rayvis
from matplotlib import pyplot as plt
import sys
import os
import msgpack
import numpy as np


def density(x):
    return 1.7145e+24 / 2 * (1 + x ** 2)


def analytic(y, y_displacement):
    eff_crit_dist = np.sqrt(2) / 2
    crit_dens = 1.7145e+24
    eff_crit_dens = density(eff_crit_dist)
    initial_density = 1.7145e+24 / 2
    omega = np.sqrt((crit_dens - initial_density) / (crit_dens - eff_crit_dens))
    print(np.pi / 2 / omega)
    return eff_crit_dist * np.sin(omega * (y - y_displacement))


def read_msgpack_energies(open_file):
    arrays = msgpack.unpackb(open_file.read())
    return [np.asarray(arr) for arr in arrays]


def generate_analytic_rays(rays, mesh):
    for ray in rays:
        y = np.linspace(ray.y[0], max(mesh.nodes[:, 1]))
        yield rayvis.Ray(analytic(y, ray.y[0]), y)


def plot_rays(folder):
    fig, axes = plt.subplots()
    with open(os.path.join(folder, "output/mesh.mfem")) as f:
        mesh = rayvis.read_mfem_mesh(f)
    with open(os.path.join(folder, "output/absorbed_energy.gf")) as f:
        grid_function = rayvis.read_grid_function(f, mesh)
    with open(os.path.join(folder, "output/rays.msgpack"), "rb") as f:
        rays = rayvis.read_msgpack_rays(f)

    poly_collection = rayvis.plot_grid_function(axes, grid_function, cmap="GnBu")
    fig.colorbar(poly_collection)
    rayvis.plot_mesh(axes, mesh, linewidth=0.5)
    rayvis.plot_rays(axes, rays, linewidth=0.5, color="red")
    analytic_rays = list(generate_analytic_rays(rays, mesh))

    plt.savefig(os.path.join(folder, "output/rays.png"))

    with open(os.path.join(folder, "output/energies0.msgpack"), "rb") as f:
        energies = read_msgpack_energies(f)

    final_energies = [e[-1] for e in energies]
    final_x = np.array([ray.x[-1] for ray in rays])
    final_analytic_x = np.array([ray.x[-1] for ray in analytic_rays])

    fig, axes = plt.subplots()
    axes.plot(final_x, final_energies, "o")
    axes.plot(final_analytic_x, final_energies, "x")

    fig, axes = plt.subplots()
    axes.plot(final_analytic_x, np.absolute(final_x - final_analytic_x))

    plt.show()


if __name__ == '__main__':
    plot_rays(sys.argv[1])

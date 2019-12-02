import json
import sys
import matplotlib.pyplot as plt
import scipy.interpolate
import matplotlib.collections
import numpy as np


def plot_quad_mesh(points, quads_as_indexes, axis=None, **kwargs):

    if not axis:
        axis = plt.gca()

    vertices = points[quads_as_indexes]
    collection = matplotlib.collections.PolyCollection(vertices, **kwargs)
    axis.add_collection(collection)
    axis.autoscale()


def plot_mesh(points, elements):

    plt.figure()
    plt.gca().set_aspect('equal')

    plot_quad_mesh(points, np.asarray(elements), color="crimson", facecolor="None")

    #plt.plot(*zip(*points), marker="o", ls="", color="crimson")

    plt.xlabel('x')
    plt.ylabel('y')


if __name__ == "__main__":
    mesh_filename = sys.argv[1]
    with open(mesh_filename) as mesh_file:
        mesh = json.load(mesh_file)
    nodes = np.asarray(mesh["points"])
    elements = np.asarray(mesh["quadrilaterals"])

    plot_mesh(nodes, elements)

    ray_filename = sys.argv[2]
    with open(sys.argv[2]) as ray_file:
        data = json.load(ray_file)

    for ray in data["rays"]:
        x = np.asarray(ray)[:, 0]
        y = np.asarray(ray)[:, 1]

    absorbed_energy = np.asarray(data["absorbedEnergy"])
    x = absorbed_energy[:, 0]
    y = absorbed_energy[:, 1]
    z = absorbed_energy[:, 2]

    print(sum(z))

    # Set up a regular grid of interpolation points
    xi, yi = np.linspace(x.min(), x.max(), 100), np.linspace(y.min(), y.max(), 100)
    xi, yi = np.meshgrid(xi, yi)

    # Interpolate
    rbf = scipy.interpolate.Rbf(x, y, z, function='linear')
    zi = rbf(xi, yi)

    #plt.imshow(zi, vmin=z.min(), vmax=z.max(), origin='lower',
    #           extent=[x.min(), x.max(), y.min(), y.max()])
    plt.scatter(x, y, c=z)
    plt.colorbar()
    plt.show()




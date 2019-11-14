import json
import sys
import matplotlib.pyplot as plt
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
        ray = json.load(ray_file)
    x = np.asarray(ray["intersections"])[:, 0]
    y = np.asarray(ray["intersections"])[:, 1]
    plt.plot(x, y, "o-")

    plt.show()






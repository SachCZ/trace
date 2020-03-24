import matplotlib.pyplot as plt
import matplotlib.collections
import numpy as np
import vtk
import json


def read_rays(filename):
    with open(filename) as rays_file:
        rays_root = json.load(rays_file)
    return rays_root["rays"]


def plot_rays(axis, rays):
    for ray in rays:
        ray = np.asarray(ray)
        x = ray[:, 0]
        y = ray[:, 1]
        axis.plot(x, y, "o-")


def read_vtk(filename):
    reader = vtk.vtkGenericDataObjectReader()
    reader.SetFileName(filename)
    reader.Update()

    nodes = np.array(reader.GetOutput().GetPoints().GetData())
    nodes = nodes[:, :-1]
    elements = np.array(reader.GetOutput().GetCells().GetData()).reshape((-1, 4))
    elements = elements[:, 1:]
    return nodes, elements


def read_grid_function(filename):
    return np.genfromtxt(filename, skip_header=5)


def plot_vtk_mesh(axis, nodes, element_indexes):
    poly_collection = matplotlib.collections.PolyCollection(nodes[element_indexes], edgecolor="blue", facecolors="")
    axis.add_collection(poly_collection)
    axis.autoscale()


def plot_grid_function(fig, axis, nodes, element_indexes, values):
    poly_collection = matplotlib.collections.PolyCollection(nodes[element_indexes], edgecolor="", cmap="Blues")
    poly_collection.set_array(values)
    fig.colorbar(poly_collection, ax=axis)
    axis.add_collection(poly_collection)
    axis.autoscale()

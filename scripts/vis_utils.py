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

    points_count = reader.GetOutput().GetCell(0).GetNumberOfPoints()

    elements = np.array(reader.GetOutput().GetCells().GetData()).reshape((-1, points_count + 1))
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


def read_gradient(filename):
    data = np.genfromtxt(filename, delimiter=",")
    x = data[:, 0]
    y = data[:, 1]
    grad_x = data[:, 2]
    grad_y = data[:, 3]
    return x, y, grad_x, grad_y


def plot_quiver(axis, x, y, f_x, f_y, **kwargs):
    u, v = np.meshgrid(f_x, f_y)
    axis.quiver(x, y, u, v, **kwargs)

def sample_analytic_gradient(x, y, function):
    grad_x = []
    grad_y = []
    for x0, y0 in zip(x,y):
        grad_x0, grad_y0 = function(x0, y0)
        grad_x.append(grad_x0)
        grad_y.append(grad_y0)
    return np.array(grad_x), np.array(grad_y)

import os
import sys

import numpy as np
import rayvis
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable


def calc_max_norm(field, analytic):
    error = field - analytic
    return max(error.norm())


def load_vector_fields(file_names):
    for file_name in file_names:
        with open(file_name, "rb") as f:
            yield rayvis.read_vector_field(f)


def calc_norms(vector_fields, analytic_vector_fields):
    for (vector_field, analytic_vector_field) in zip(vector_fields, analytic_vector_fields):
        yield calc_max_norm(vector_field, analytic_vector_field)


def plot_grad_summary(fig, dual_mesh, grid_function, gradient, analytic):
    gs = GridSpec(nrows=2, ncols=2, width_ratios=[1, 1], height_ratios=[1, 1])
    fun_axes = fig.add_subplot(gs[0, 0])
    grad_axes = fig.add_subplot(gs[1, 0])
    analytic_axes = fig.add_subplot(gs[0, 1])
    diff_axes = fig.add_subplot(gs[1, 1])

    contour = rayvis.plot_vector_field(grad_axes, gradient, dual_mesh)
    fig.colorbar(contour, cax=make_axes_locatable(grad_axes).append_axes("right", size="5%", pad=0.05))

    contour = rayvis.plot_grid_function(fun_axes, grid_function)
    fig.colorbar(contour, cax=make_axes_locatable(fun_axes).append_axes("right", size="5%", pad=0.05))

    error = gradient - analytic
    grad_norm = np.asarray(2 * [gradient.norm()])
    contour = rayvis.plot_vector_field(diff_axes, error, dual_mesh)
    fig.colorbar(contour, cax=make_axes_locatable(diff_axes).append_axes("right", size="5%", pad=0.05))

    contour = rayvis.plot_vector_field(analytic_axes, analytic, dual_mesh)
    analytic_axes.set_xlabel("analytic")
    fig.colorbar(contour, cax=make_axes_locatable(analytic_axes).append_axes("right", size="5%", pad=0.05))


def main(folder):
    path = os.path.join(folder, "output/grad{}.msgpack")
    analytic_path = os.path.join(folder, "output/analytic_grad{}.msgpack")
    segments_count = range(10, 21, 10)
    file_names = [path.format(count) for count in segments_count]
    analytic_file_names = [analytic_path.format(count) for count in segments_count]
    vector_fields = list(load_vector_fields(file_names))
    analytic_vector_fields = list(load_vector_fields(analytic_file_names))
    norms = calc_norms(vector_fields, analytic_vector_fields)
    fig, axes = plt.subplots()
    h_eff = np.asarray(1) / segments_count
    axes.loglog(h_eff, list(norms), "o")
    axes.loglog(h_eff, h_eff ** 2)
    plt.savefig(os.path.join(folder, "output/conv.png"))

    for segments, vector_field, analytic_vector_field in zip(segments_count, vector_fields, analytic_vector_fields):
        fig = plt.figure(figsize=(15, 7.5), dpi=300)
        with open(os.path.join(folder, "output/mesh{}.mfem".format(segments))) as f:
            mesh = rayvis.read_mfem_mesh(f)
        with open(os.path.join(folder, "output/dual_mesh{}.mfem".format(segments))) as f:
            dual_mesh = rayvis.read_mfem_mesh(f)
        with open(os.path.join(folder, "output/func{}.gf".format(segments))) as f:
            gf = rayvis.read_grid_function(f, mesh)
        plot_grad_summary(fig, dual_mesh, gf, vector_field, analytic_vector_field)
        plt.savefig(os.path.join(folder, "output/summary{}.png".format(segments)))


if __name__ == '__main__':
    main(sys.argv[1])

import os
import sys

import numpy as np
import rayvis
from matplotlib import pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.axes_grid1 import make_axes_locatable
import yaml
import scipy.optimize as optimize


def calc_max_norm(field, analytic, h):
    error = field - analytic
    return max(error.norm()) / max(analytic.norm())


def calc_l2_norm(field, analytic, h):
    error = field - analytic
    return np.sqrt(np.sum(error.norm() ** 2)) / np.sqrt(np.sum(analytic.norm() ** 2))


def load_vector_fields(file_names):
    for file_name in file_names:
        with open(file_name, "rb") as f:
            yield rayvis.read_vector_field(f)


def calc_norms(vector_fields, analytic_vector_fields, norm_func):
    for (vector_field, analytic_vector_field) in zip(vector_fields, analytic_vector_fields):
        yield norm_func(vector_field, analytic_vector_field, 1.0 / (len(vector_field.coord_x) - 1))


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


def load_config(filename):
    with open(filename, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def main(folder, config_name="config.yaml"):
    path = os.path.join(folder, "output/grad{}_{}.msgpack")
    analytic_path = os.path.join(folder, "output/analytic_grad{}_{}.msgpack")
    config = load_config(os.path.join(folder, "input", config_name))
    segments_count = [int(count) for count in config["meshes"]["segments"]]
    factors = config["meshes"]["random_factors"]
    conv_fig, conv_axes = plt.subplots()

    displacement_from0 = []

    for factor in factors:
        file_names = [path.format(count, factor) for count in segments_count]
        analytic_file_names = [analytic_path.format(count, factor) for count in segments_count]
        vector_fields = list(load_vector_fields(file_names))
        analytic_vector_fields = list(load_vector_fields(analytic_file_names))
        norms = list(calc_norms(vector_fields, analytic_vector_fields, calc_l2_norm))
        h_eff = np.asarray(1) / segments_count

        try:
            def fit_func(h, e, a, b):
                return a * h ** e + b

            popt, pcov = optimize.curve_fit(fit_func, h_eff, norms, p0=[1., 1., 0])
            b = popt[2]
            a = popt[1]
            exponent = popt[0]
        except RuntimeError as e:
            def fit_func(h, a):
                return a * np.ones(len(h))

            popt, pcov = optimize.curve_fit(fit_func, h_eff, norms, p0=[1.])
            b = popt[0]
            a = 0
            exponent = 0

        line, = conv_axes.plot(h_eff, norms, "o", label="factor = {}, $f(h) = {:.0f}x^{{{:.2f}}} + {:.3f}$".format(factor, a, exponent, b))
        x = np.linspace(min(h_eff), max(h_eff), 100)
        conv_axes.plot(x, fit_func(x, *popt), color=line.get_color())
        displacement_from0.append(b)

        for segments, vector_field, analytic_vector_field in zip(segments_count, vector_fields, analytic_vector_fields):
            if segments > 30:
                continue
            fig = plt.figure(figsize=(10, 5), dpi=150)
            with open(os.path.join(folder, "output/mesh{}_{}.mfem".format(segments, factor))) as f:
                mesh = rayvis.read_mfem_mesh(f)
            with open(os.path.join(folder, "output/dual_mesh{}_{}.mfem".format(segments, factor))) as f:
                dual_mesh = rayvis.read_mfem_mesh(f)
            with open(os.path.join(folder, "output/func{}_{}.gf".format(segments, factor))) as f:
                gf = rayvis.read_grid_function(f, mesh)
            plot_grad_summary(fig, dual_mesh, gf, vector_field, analytic_vector_field)
            plt.savefig(os.path.join(folder, "output/summary{}_{}.png".format(segments, factor)))
    conv_axes.legend()
    conv_fig.savefig(os.path.join(folder, "output/conv.png"))

    conv_error_fig, conv_error_axes = plt.subplots()
    conv_error_axes.plot(factors, displacement_from0, "o")
    conv_error_axes.grid()
    conv_error_axes.set_xlabel("factor")
    conv_error_axes.set_ylabel("$b$")
    conv_error_fig.savefig(os.path.join(folder, "output/conv_error.png"))


if __name__ == '__main__':
    args = list(sys.argv)
    if len(args) == 2:
        main(sys.argv[1])
    elif len(sys.argv) == 3:
        main(sys.argv[1], sys.argv[2])
    else:
        print("Invalid args")

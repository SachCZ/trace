import numpy as np
import rayvis
from matplotlib import pyplot as plt
import yaml
import scipy.optimize as optimize


def calc_max_norm(field, analytic):
    error = field - analytic
    return max(error.norm()) / max(analytic.norm())


def calc_l2_norm(field, analytic):
    error = field - analytic
    return np.sqrt(np.sum(error.norm() ** 2)) / np.sqrt(np.sum(analytic.norm() ** 2))


def load_vector_fields(file_names):
    for file_name in file_names:
        with open(file_name, "rb") as f:
            yield rayvis.read_vector_field(f)


def calc_norms(vector_fields, analytic_vector_fields, norm_func):
    for (vector_field, analytic_vector_field) in zip(vector_fields, analytic_vector_fields):
        yield norm_func(vector_field, analytic_vector_field)


def load_config(filename):
    with open(filename, 'r') as stream:
        try:
            return yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            print(exc)


def fit_func(h, a):
    return a * h ** -2


def factors_to_norms(factors, segments_counts):
    for factor in factors:
        file_names = ["output/grad{}_{}.msgpack".format(count, factor) for count in segments_counts]
        analytic_file_names = ["output/analytic_grad{}_{}.msgpack".format(count, factor) for count in segments_counts]
        vector_fields = list(load_vector_fields(file_names))
        analytic_vector_fields = list(load_vector_fields(analytic_file_names))
        yield list(calc_norms(vector_fields, analytic_vector_fields, calc_l2_norm))


def main():
    segments_counts = [8, 16, 32, 64, 128, 256, 512]
    factors = [0, 0.02, 0.04]
    markers = ["x", "+", "1"]
    conv_fig, conv_axes = plt.subplots()
    norms_set = list(factors_to_norms(factors, segments_counts))

    popt, pcov = optimize.curve_fit(fit_func, segments_counts, norms_set[0], p0=[14.])
    x = np.linspace(min(segments_counts), max(segments_counts), 100)
    conv_axes.loglog(x, fit_func(x, *popt), label="$~x^{-2}$")

    for factor, norms, marker in zip(factors, norms_set, markers):
        conv_axes.loglog(segments_counts, norms, marker, label="$f$ = {}".format(factor))

    conv_axes.legend()
    conv_axes.grid()
    conv_axes.set_xlabel("$N$ [-]")
    conv_axes.set_ylabel("$\\Delta \\varepsilon / \\varepsilon$ [-]")
    plt.show()
    conv_fig.savefig("images/conv.png")

    return norms_set


if __name__ == '__main__':
    main()

import numpy as np
import rayvis
from matplotlib import pyplot as plt
import yaml
import scipy.optimize as optimize


def calc_l2_norm(field, analytic):
    error = field - analytic
    return np.sqrt(np.sum(error.norm() ** 2)) / np.sqrt(np.sum(analytic.norm() ** 2))


def load_vector_fields(file_names):
    for file_name in file_names:
        with open(file_name, "rb") as f:
            yield rayvis.read_vector_field(f)


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
        yield [calc_l2_norm(vf, avf) for (vf, avf) in zip(vector_fields, analytic_vector_fields)]


def main():

    config = load_config("input/config.yaml")
    segments_counts = config["meshes"]["segments"]
    factors = config["meshes"]["random_factors"]
    markers = ["x", "+", "1"]
    conv_fig, conv_axes = plt.subplots()
    norms_set = list(factors_to_norms(factors, segments_counts))

    for factor, norms, marker in zip(factors, norms_set, markers):
        conv_axes.semilogx(segments_counts, norms, marker, label="$f$ = {}".format(factor))



    conv_axes.grid()
    conv_axes.set_xlabel("$N$ [-]")
    conv_axes.set_ylabel("$\\Delta \\varepsilon / \\varepsilon$ [-]")

    plt.legend(bbox_to_anchor=(0.1, 0.6), bbox_transform=conv_axes.transAxes)

    plt.show()
    conv_fig.savefig("images/green_lin_conv.png")


if __name__ == '__main__':
    main()

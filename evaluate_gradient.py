import numpy as np
import rayvis
from matplotlib import pyplot as plt


def analytic_fun(x, y):
    return 3.69106 * np.cos(11.9066 * x + 9.36195 * y), 2.9022 * np.cos(11.9066 * x + 9.36195 * y)


def calc_max_norm(field, analytic):
    error = field - analytic(field.coord_x, field.coord_y)
    return max(error.norm())


def load_vector_fields(file_names):
    for file_name in file_names:
        with open(file_name, "rb") as f:
            yield rayvis.read_vector_field(f)


def calc_norms(vector_fields):
    for vector_field in vector_fields:
        yield calc_max_norm(vector_field, analytic_fun)


def main():
    path = "run/gradient/output/household{}.msgpack"
    segments_count = range(10, 101, 10)
    file_names = [path.format(count) for count in segments_count]
    vector_fields = list(load_vector_fields(file_names))
    norms = calc_norms(vector_fields)
    fig, axes = plt.subplots()
    h_eff = np.asarray(1) / segments_count
    axes.loglog(h_eff, list(norms), "o")
    axes.loglog(h_eff, h_eff ** 2)

    fig, axes = plt.subplots()
    rayvis.plot_vector_field(axes, vector_fields[1])
    plt.show()


if __name__ == '__main__':
    main()

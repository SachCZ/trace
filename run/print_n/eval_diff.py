from copy import deepcopy

import rayvis
import matplotlib as mpl
from matplotlib import pyplot as plt
import scipy.integrate as integrate
import numpy as np
from matplotlib.patches import Rectangle


def linear_fun(n, n0, n1, z0, z1):
    k = (n0 - n1) / (z1 - z0)
    b = n1 + k * z1
    return (n - b) / k


def solution(z0, zw, n0, nw, nc):
    zc = linear_fun(nc, nw, n0, zw, z0)
    return n0 * (zc - zw) * nc * (1 - nw ** 2)


def trajectory(y, x0, y0, grad, ne, nc):
    return x0 - grad / 4 / (nc - ne) * (y - y0) ** 2


def crit_dens(lamb):
    me = 9.1093897e-28
    c = 2.99792458e10
    e = 4.8032068e-10
    return me * np.pi * c ** 2 / e ** 2 / lamb ** 2


if __name__ == '__main__':
    with open("input/mesh.vtk") as f:
        mesh = rayvis.read_vtk_mesh(f)
        mesh.elements = mesh.elements[mesh.elements[:, 0].argsort()]

    with open("input/gain.gf") as f:
        gain = rayvis.read_grid_function(f, mesh)

    x_orig = np.asarray(
        [0.02216, 0.0156218, 0.0124454, 0.0100609, 0.00811045, 0.0064038, 0.00478522, 0.00346796, 0.00240774, 0.0015729,
         0.00094264, 0.000514716, 0.00026871, 0.000139284, 7.17372e-05, 3.61921e-05, 9.20324e-06, -1.17732e-05,
         -2.97418e-05, -4.49682e-05, -5.82379e-05, -6.9885e-05, -8.02068e-05, -8.93425e-05, -9.73705e-05, -0.000104184,
         -0.000109724, -0.000115328, -0.000120954, -0.000126583, -0.000132228, -0.000137906, -0.000143703, -0.000149669,
         -0.000156099, -0.000163073, -0.000170687, -0.000179149, -0.000188746, -0.000199971, -0.000214085, -0.000233533,
         -0.000259081, -0.000287604, -0.000319227, -0.000354285, -0.000393152, -0.000436241, -0.000484013, -0.000536974,
         -0.00059569, -0.000660785, -0.000732953, -0.000812961, -0.000901662, -0.001]) * 1e4
    x = deepcopy(x_orig)
    size_y = len(gain.values) / (len(x) - 1)
    gain = gain.values.reshape((int(size_y), len(x) - 1))
    gain = gain[59, :]
    ne = np.asarray(
        [2.98409e+19, 5.23045e+19, 9.84547e+19, 1.44417e+20, 1.93469e+20, 2.38541e+20, 3.07437e+20, 4.23456e+20,
         5.87307e+20, 8.30304e+20, 1.22976e+21, 1.8705e+21, 1.65413e+21, 1.36323e+21, 3.79627e+21, 3.62981e+21,
         5.43067e+21, 6.46437e+21, 7.51983e+21, 9.0897e+21, 1.09211e+22, 1.34127e+22, 1.69354e+22, 2.23067e+22,
         3.14185e+22, 4.81108e+22, 6.40848e+22, 7.39701e+22, 8.60965e+22, 1.00358e+23, 1.16529e+23, 1.33343e+23,
         1.49995e+23, 1.61904e+23, 1.68126e+23, 1.72918e+23, 1.74557e+23, 1.71239e+23, 1.61665e+23, 1.41441e+23,
         1.09439e+23, 8.1022e+22, 7.03497e+22, 6.9965e+22, 6.9965e+22, 6.99649e+22, 6.99649e+22, 6.99649e+22,
         6.99649e+22, 6.99649e+22, 6.99649e+22, 6.99649e+22, 6.99649e+22, 6.99649e+22, 6.99649e+22])
    n = np.asarray(
        [0.999991, 0.999985, 0.999971, 0.999958, 0.999944, 0.99993, 0.99991, 0.999876, 0.999829, 0.999758, 0.999641,
         0.999454, 0.999517, 0.999602, 0.998893, 0.998942, 0.998419, 0.998121, 0.99782, 0.997375, 0.996868, 0.996198,
         0.995292, 0.994003, 0.992033, 0.988949, 0.986547, 0.985328, 0.983977, 0.982547, 0.981061, 0.979629, 0.978313,
         0.977479, 0.97714, 0.976931, 0.976966, 0.977368, 0.978248, 0.98004, 0.983131, 0.986461, 0.988021, 0.988163,
         0.988163, 0.988165, 0.988166, 0.988166, 0.988166, 0.988166, 0.988166, 0.988166, 0.988166, 0.988166, 0.988166])
    lamb = 25.5e-7
    nc = crit_dens(lamb)
    delta_x = np.absolute(x[1:] - x[:-1])
    x = x[1:]
    grad = (ne[1:] - ne[:-1]) / (x[1:] - x[:-1])

    fig, [ax0, ax_zoom] = plt.subplots(1, 2)
    xx = []
    ne_xx = []
    y_maxes = []
    for x_i, dx, grad_i, ne_i in zip(x[:-1], delta_x, grad, ne[:-1]):
        if x_i < 2:
            continue
        if grad_i >= 0:
            continue
        y_max = np.sqrt(4 * (nc - ne_i) / -grad_i * dx)
        y_maxes.append(y_max)
        xx.append(x_i)
        ne_xx.append(ne_i)
        y = np.linspace(0, y_max)
        print(x_i, dx, y_max, "dep: ", dx / y_max)
        ax0.add_patch(Rectangle((5, 0), 44, 2900, facecolor="none", edgecolor="red", linewidth=1))

        ax0.plot(trajectory(y, x_i, 0, grad_i, ne_i, nc), y, color="white")
        ax_zoom.plot(trajectory(y, x_i, 0, grad_i, ne_i, nc), y, color="white")

    X, Y = np.meshgrid([-10, 26000, 27000], [x_orig[0]] + list(xx))
    Z = np.asarray(list(zip(ne_xx, ne_xx)))
    cmap = ax0.pcolormesh(Y, X, Z, shading='flat', edgecolor="black")
    plt.colorbar(cmap, ax=ax_zoom)
    ax_zoom.pcolormesh(Y, X, Z, shading='flat', edgecolor="black")
    ax0.set_ylim((0, 25000))
    ax0.set_xlim((5, x_orig[0]))
    ax_zoom.set_ylim((0, 2900))
    ax_zoom.set_xlim((5, 47))
    plt.setp(ax_zoom.spines.values(), color="red", linewidth=2)
    plt.tight_layout()
    ax_zoom.set_xlabel("x[μm]")
    ax_zoom.set_ylabel("y[μm]")
    ax0.set_xlabel("x[μm]")
    ax0.set_ylabel("y[μm]")
    ax0.set_title("$n_\\mathrm{e}$ [cm$^{-3}$]")
    ax_zoom.set_title("$n_\\mathrm{e}$ [cm$^{-3}$]")
    plt.tight_layout()
    plt.savefig("diffr2d.png")
    plt.show()

    fig, host = plt.subplots()
    ax2 = host.twinx()
    p1, = host.plot(xx, y_maxes, color="red", marker="+", label="diffraction scale")
    host.grid()
    p2, = ax2.plot(x, gain[::-1], color="blue", linestyle="--", marker="x", label="$gain coefficient")
    host.yaxis.label.set_color(p1.get_color())
    ax2.yaxis.label.set_color(p2.get_color())
    host.tick_params(axis='y', colors=p1.get_color())
    ax2.tick_params(axis='y', colors=p2.get_color())
    host.set_xlabel("x[μm]")
    host.set_ylabel("y[μm]")
    ax2.set_ylabel("$g$ [cm$^{-1}$]")
    lns = [p1, p2]
    host.legend(handles=lns, loc=9)
    plt.savefig("diffr_gain.png")
    plt.show()

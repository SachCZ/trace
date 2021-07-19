import rayvis
from matplotlib import pyplot as plt
import numpy as np


def density(x):
    return 1.7145e+24 / 2 * (1 + x ** 2)


def analytic(y, y_displacement):
    eff_crit_dist = np.sqrt(2) / 2
    crit_dens = 1.7145e+24
    eff_crit_dens = density(eff_crit_dist)
    initial_density = 1.7145e+24 / 2
    omega = np.sqrt((crit_dens - initial_density) / (crit_dens - eff_crit_dens))
    return eff_crit_dist * np.sin(omega * (y - y_displacement))


def compare_errors(axes, grad_type, **kwargs):
    with open("output/rays512{}0.02.msgpack".format(grad_type), "rb") as f:
        rays = rayvis.read_msgpack_rays(f)

    final_x = np.array([ray.x[-1] for ray in rays])
    final_analytic_x = np.array([analytic(ray.y[-1], ray.y[0]) for ray in rays])

    error = np.absolute(final_x - final_analytic_x)

    axes.plot(final_x, error, linestyle="", **kwargs)


def main():
    fig, [ax0, ax1] = plt.subplots(2)
    compare_errors(ax0, "ls", marker="+", label="ls")
    compare_errors(ax0, "mfem", marker=".", label="mfem")
    compare_errors(ax0, "integral", marker="1", label="integral")
    ax0.set_xlabel("$x$ [cm]")
    ax0.set_ylabel("$|x - x_{ref}|$ [cm]")
    ax0.set_ylim((0, 0.0008))
    ax0.legend()
    ax0.grid()

    with open("output/rays512ls0.02.msgpack", "rb") as f:
        rays_ls = rayvis.read_msgpack_rays(f)

    rayvis.plot_rays(ax1, rays_ls, linewidth=0.5, color="red")

    ax1.set_xlabel("$x$ [cm]")
    ax1.set_ylabel("$y$ [cm]")
    ax1.grid()

    plt.tight_layout()
    plt.savefig("images/error512.png")


if __name__ == '__main__':
    main()

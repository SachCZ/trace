import rayvis
from matplotlib import pyplot as plt
import numpy as np
import scipy.optimize as optimize


def density(x):
    return 1.7145e+24 / 2 * (1 + x ** 2)


def analytic(y, y_displacement):
    eff_crit_dist = np.sqrt(2) / 2
    crit_dens = 1.7145e+24
    eff_crit_dens = density(eff_crit_dist)
    initial_density = 1.7145e+24 / 2
    omega = np.sqrt((crit_dens - initial_density) / (crit_dens - eff_crit_dens))
    return eff_crit_dist * np.sin(omega * (y - y_displacement))


def fit_func(x, a, b):
    return a * x ** b


def compare_gradients(axes, first_grad_type, second_grad_type, eval_method=np.median, fit=False):
    labels = {
        "ls": "LSQ",
        "mfem": "FEM",
        "integral": "Green's formula"
    }

    for grad_type, markers in zip([first_grad_type, second_grad_type], [["x", "+", "."], ["1", "2", "3"]]):
        for factor, mark in zip([0, 0.02, 0.04], markers):
            if second_grad_type == "mfem" and grad_type == "ls" and factor > 0:
                continue

            errors = []
            counts = [8, 16, 32, 64, 128, 256, 512]
            for count in counts:
                with open("output/rays{}{}{}.msgpack".format(count, grad_type, factor), "rb") as f:
                    rays = rayvis.read_msgpack_rays(f)

                final_x = np.array([ray.x[-1] for ray in rays])
                final_analytic_x = np.array([analytic(ray.y[-1], ray.y[0]) for ray in rays])

                error = np.absolute(final_x - final_analytic_x)
                errors.append(eval_method(error))

            axes.loglog(counts, errors, mark, label=labels[grad_type] + " $f = $" + str(factor))

            if fit and grad_type == "ls" and factor == 0:
                popt, pcov = optimize.curve_fit(fit_func, counts, errors, p0=[1, -2])
                x = np.linspace(min(counts), max(counts))
                axes.loglog(x, fit_func(x, *popt), label="$\\sim N^{-1.87}$", color="blue")
                print(popt[1])

    handles, labels = axes.get_legend_handles_labels()

    handles = [handles[1], handles[0], *handles[2:]]
    labels = [labels[1], labels[0], *labels[2:]]

    axes.legend(handles, labels, loc=3)
    axes.grid()


def main():
    fig, [ax0, ax1] = plt.subplots(1, 2, sharey="row", figsize=(10.0, 4.8))
    compare_gradients(ax0, "ls", "integral", fit=True)
    ax0.set_xlabel("$N$ [-]")
    ax0.set_ylabel("median($|x - x_{ref}|$) [cm]")
    compare_gradients(ax1, "ls", "mfem", fit=True)
    ax1.set_xlabel("$N$ [-]")
    plt.tight_layout()
    plt.savefig("images/error_median.png")
    plt.close()

    fig, [ax0, ax1] = plt.subplots(1, 2, sharey="row", figsize=(10.0, 4.8))
    compare_gradients(ax0, "ls", "integral", np.max)
    ax0.set_xlabel("$N$ [-]")
    ax0.set_ylabel("max($|x - x_{ref}|$) [cm]")
    compare_gradients(ax1, "ls", "mfem", np.max)
    ax1.set_xlabel("$N$ [-]")
    plt.tight_layout()
    plt.savefig("images/error_max.png")
    plt.close()


if __name__ == '__main__':
    main()

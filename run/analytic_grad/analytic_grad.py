import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    fig, [u_ax, grad_u_ax] = plt.subplots(2)
    x = np.linspace(0, 1, 1000)
    y = np.linspace(0, 1, 1000)
    x, y = np.meshgrid(x, y)
    u = 1.62 + 0.31 * np.sin(np.pi * (3.79 * x + 2.98 * y))
    grad_u = (1.1749 * np.pi * np.cos(np.pi * (3.79 * x + 2.98 * y)), 0.9238 * np.pi * np.cos(np.pi * (3.79 * x + 2.98 * y)))
    grad_u_norm = np.sqrt(grad_u[0]**2 + grad_u[1]**2)
    cmap = u_ax.pcolormesh(x, y, u, shading="auto", cmap="plasma")
    plt.colorbar(cmap, ax=u_ax)
    u_ax.set_xlabel("x")
    u_ax.set_ylabel("y")
    cmap = grad_u_ax.pcolormesh(x, y, grad_u_norm)
    plt.colorbar(cmap, ax=grad_u_ax)
    grad_u_ax.set_xlabel("x")
    grad_u_ax.set_ylabel("y")
    plt.savefig("images/analytic.png")
    plt.show()

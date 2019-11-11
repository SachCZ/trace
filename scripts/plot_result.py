import json
import sys

import numpy as np
from matplotlib import pyplot as plt

if __name__ == "__main__":

    mesh_filename = sys.argv[1]
    with open(mesh_filename) as mesh_file:
        mesh = json.load(mesh_file)
    x = np.asarray(mesh["points"])[:, 0]
    y = np.asarray(mesh["points"])[:, 1]
    triangles = np.asarray(mesh["triangles"])
    plt.triplot(x, y, triangles, '-')

    ray_filename = sys.argv[2]
    with open(sys.argv[2]) as ray_file:
        ray = json.load(ray_file)
    x = np.asarray(ray["intersections"])[:, 0]
    y = np.asarray(ray["intersections"])[:, 1]
    plt.plot(x, y, "o-")

    plt.show()




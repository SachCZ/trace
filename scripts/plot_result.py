#!/usr/bin/env python3

import json
import sys
import matplotlib.pyplot as plt
import numpy as np

if __name__ == "__main__":
    rays_filename = sys.argv[1]
    with open(rays_filename) as rays_file:
        rays_root = json.load(rays_file)
    rays = rays_root["rays"]
    for ray in rays:
        ray = np.asarray(ray)
        x = ray[:, 0]
        y = ray[:, 1]

        plt.plot(x, y, "o-")

    plt.axis('equal')
    plt.show()

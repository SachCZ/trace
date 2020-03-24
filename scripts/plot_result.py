#!/usr/bin/env python3

import sys
import matplotlib.pyplot as plt
import scripts.vis_utils as vis_utils

if __name__ == "__main__":
    _, axis = plt.subplots()
    vis_utils.plot_rays(axis, vis_utils.read_rays(sys.argv[1]))

    nodes, elements = vis_utils.read_vtk(sys.argv[2])
    vis_utils.plot_vtk_mesh(axis, nodes, elements)

    plt.axis('equal')
    plt.show()

    fig, axis = plt.subplots()
    values = vis_utils.read_grid_function(sys.argv[3])
    vis_utils.plot_grid_function(fig, axis, nodes, elements, values)
    #vis_utils.plot_rays(axis, vis_utils.read_rays(sys.argv[1]))
    plt.axis('equal')
    plt.show()

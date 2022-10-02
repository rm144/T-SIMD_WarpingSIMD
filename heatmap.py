# ===========================================================================
# 
# heatmap.py --
# Christoph Berganski
# 
# This source code file is part of the following software:
# 
#    - the low-level C++ template SIMD library
#    - the SIMD implementation of the MinWarping and the 2D-Warping methods 
#      for local visual homing.
# 
# The software is provided based on the accompanying license agreement
# in the file LICENSE or LICENSE.doc. The software is provided "as is"
# without any warranty by the licensor and without any liability of the
# licensor, and the software may not be distributed by the licensee; see
# the license agreement for details.
# 
# (C) Ralf Möller
#     Computer Engineering
#     Faculty of Technology
#     Bielefeld University
#     www.ti.uni-bielefeld.de
# 
# ===========================================================================

# Plotting, saving figures and images, ...
import matplotlib.pyplot as plt
# Arrays and images
import numpy as np


# This is to be able to edit text in SVG files later
# plt.rcParams['svg.fonttype'] = 'none'


# Creates a color mesh heatmap plot
def heatmap(data: np.array, title: str = 'Heatmap', argmin=None, estmin=False):
    # Decompose columns of the search log
    x, y, alpha, psi, d = data.T
    # Width and height of the x-y range
    width, height = len(np.unique(x)), len(np.unique(y))
    # Reshape to recover the 2d search range
    x = x.reshape((height, width))
    y = y.reshape((height, width))
    d = d.reshape((height, width))
    # Plot a grayscale heatmap of the warping d values over the (x,y) search
    # range
    plt.pcolormesh(x, y, d, cmap='gray')
    # Mark the true minimum
    if argmin is not None:
        plt.scatter(argmin[0], argmin[1], color='green', s=500, marker='+')
    # Plot a blue marker at the search minimum
    if estmin:
        # Find the location of the estimated minimum
        estmin = [
            x[np.unravel_index(d.argmin(), (height, width))],
            y[np.unravel_index(d.argmin(), (height, width))]
        ]
        plt.scatter(*estmin, edgecolors='red', facecolors='blue', s=75)
    # Limit the plotting range
    plt.xlim(x.min(), x.max())
    plt.ylim(y.min(), y.max())
    # Make scaling of the x- and y-axis the same
    plt.gca().set_aspect('equal')
    # Add title to plot naming the input file
    plt.title(title)
    # Add axis labels
    plt.xlabel('deltaX [rad]')
    plt.ylabel('deltaY [rad]')


# Main function to not have all the code in global scope
def main(filename: str, output: str, title: str = None, argmin=None):
    # Load the grid search log file
    log = np.loadtxt(filename)
    # Fix the figure size
    # plt.figure(figsize=(5, 5))
    # Create the heatmap plot
    heatmap(log, filename if title is None else title, argmin, estmin=True)
    # Put some space between the subplots to improve readability
    plt.tight_layout()
    # Save the plot to the output file
    plt.savefig(output)


# Entrypoint if run as script/module
if __name__ == '__main__':
    # Command line argument parsing
    import argparse

    # Set up a new argument parser
    parser = argparse.ArgumentParser()
    # Positional arguments (filename of search log to visualize)
    parser.add_argument('filename', type=str)
    parser.add_argument('output', type=str)
    parser.add_argument('--argmin', type=float, nargs=2)
    parser.add_argument('--title', type=str)
    # Parse supplied arguments using the parser
    args = parser.parse_args().__dict__

    # Run main function to not have the whole script in global scope
    main(**args)

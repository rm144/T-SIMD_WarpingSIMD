# ===========================================================================
# 
# simplex.py --
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

# Reuse heatmap plotting for background
from heatmap import heatmap


# This is to be able to edit text in SVG files later
# plt.rcParams['svg.fonttype'] = 'none'


# Main function to not have all the code in global scope
def main(filename: str, output: str, title: str = None, background: str = None,
         argmin=None):
    # Load the background if given
    background = np.loadtxt(background) if background is not None else None
    # Load the search log file
    log = np.loadtxt(filename)
    # Decompose columns of the search log
    x1s, y1s, _, _, _, x2s, y2s, _, _, _, x3s, y3s, _, _, _, = log.T
    # Rearrange coordinates to list the simplex points
    points = zip(zip(x1s, y1s), zip(x2s, y2s), zip(x3s, y3s))
    # Make the figure big enough to fit multiple plots next to each other
    plt.figure(figsize=(10, 10))
    # Add title to plot naming the input file
    plt.suptitle(filename if title is None else title)
    # Iterate over the simplex steps
    for step, simplex in enumerate(points):
        # Plot each step of the pattern into a new subplot
        plt.subplot(int(np.ceil(len(x1s) / 3)), 3, step + 1)
        # Fix and limit the plotting range
        plt.xlim(min(x1s.min(), x2s.min(), x3s.min()),
                 max(x1s.max(), x2s.max(), x3s.max()))
        plt.ylim(min(y1s.min(), y2s.min(), y3s.min()),
                 max(y1s.max(), y2s.max(), y3s.max()))
        # Plot the background if given
        if background is not None:
            heatmap(background, title=f'Step {step + 1}', argmin=argmin)
        # Plot the lines connecting the pattern points
        plt.triplot(*zip(*simplex), color='red')
        # Plot the pattern points
        plt.scatter(*zip(*simplex), color='red')
        # Make scaling of the x- and y-axis the same
        plt.gca().set_aspect('equal')
    # Plot the last best point in blue
    plt.scatter(x1s[-1], y1s[-1], edgecolors='red', facecolors='blue')
    # Put some space between the subplots to improve readability
    plt.tight_layout(pad=2.5)
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
    parser.add_argument('background', type=str, default=None, nargs='?')
    parser.add_argument('--argmin', type=float, nargs=2)
    parser.add_argument('--title', type=str)
    # Parse supplied arguments using the parser
    args = parser.parse_args().__dict__

    # Run main function to not have the whole script in global scope
    main(**args)

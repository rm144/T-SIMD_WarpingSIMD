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
    xs, ys, alpha, psi, d, wxs, wys = log.T
    # Make the figure big enough to fit multiple plots next to each other
    plt.figure(figsize=(10, 10))
    # Add title to plot naming the input file
    plt.suptitle(filename if title is None else title)
    # Iterate over the pattern steps
    for step, (x, y, wx, wy) in enumerate(zip(xs, ys, wxs, wys)):
        # Plot each step of the pattern into a new subplot
        plt.subplot(int(np.ceil(len(xs) / 3)), 3, step + 1)
        # Fix and limit the plotting range
        plt.xlim(xs.min() - wxs.max(), xs.max() + wxs.max())
        plt.ylim(ys.min() - wys.max(), ys.max() + wys.max())
        # Plot the background if given
        if background is not None:
            heatmap(background, title=f'Step {step + 1}', argmin=argmin)
        # List all pattern points +/- wx/wy along axes
        pattern = [
            (x - wx, y), (x + wx, y), (x, y - wy), (x, y + wy), (x, y)
        ]
        # Plot the lines connecting the pattern points
        plt.triplot(*zip(*pattern), color='red', mask=4 * [True])
        # Plot the pattern points
        plt.scatter(*zip(*pattern), color='red')
        # Make scaling of the x- and y-axis the same
        plt.gca().set_aspect('equal')
    # Plot the last center point (best) in blue
    plt.scatter(xs[-1], ys[-1], edgecolors='red', facecolors='blue')
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

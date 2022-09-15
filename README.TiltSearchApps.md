# Warping SIMD Tilt Correction+Search Apps

## IdealScale / idealscale

The `IdealScale` application determines the ideal scaling (pixel- and post
scale) to be used by the warping scale plane stack computation. It does this for
a single image from a panoramic image database printing the result to the
standard output:

```shell
# Usage: IdealScale root filename
$ IdealScale data/lab/d 000000
9.85177  170.664
```

Before doing any warping runs, set the ideal scale by taking the parameter-wise
maximum over all or at least a significant portion of images from the test
dataset.

idealscale is a script which determines the minimum of pixel- and postscale over an entire database for different database variants:

```shell
# Usage: idealscale database variants
$ idealscale data/lab "n d"
6.80378 170.664
```

## Tilt

The `Tilt` application reads in a panoramic image and its ground truth tilt
parameters and applies tilt correction of the selected solution and
interpolation to the image writing the resulting image to file (into the
directory specified by the third argument `output_root`, which must already
exist) and printing time measurements of the tilt correction (and the input
image descriptor) to the standard output:

```shell
# Usage: Tilt root filename output_root solution [interpolation=NEAREST]
$ Tilt data/lab/d 000000 out 0 0
000000 0 0 0.087266 -0.087266 0 7856 # Last number is runtime in microseconds
```

For comparisons, the original image is transferred to the output directory as
well with `*.orig.pgm` added to the filename.

To get meaningful runtime statistics, this should be repeated and averaged over
multiple images. Comparison of the tilt corrections can be done by applying the
different solutions to the same input image and visually inspecting the image
differences or computing an average pixel-wise difference.

## BaselineDemo

The `BaselineDemo` application demonstrates the warping homing performance on
tilted panoramic images without or with ground truth tilt correction (this means
there is no tilt search involved). The application processes a single snapshot
current view image pair with optional orientation change. All other (warping and
tilt correction related parameters) are set via environment variables. The
results are written to the standard output:

```shell
# Usage: BaselineDemo ss_root ss cv_root cv [rotate-ss] [rotate-cv] [no-tilt]
$ BaselineDemo data/lab/d 000000 data/lab/d 001000
000000 ... 001000 ... 0.05236 0.087266 0.0981748 0 40246 1 19671 22155 
# EST:   beta=3.043418, alpha=0.098175, psi=0.000000
# TRUE:  beta=3.141593, alpha=0.000000, psi=0.000000
# ERROR: beta=0.098175, alpha=0.098175, psi=0.000000
# TILT EST:   delta_X=0.052360, delta_Y=0.087266
# TILT TRUE:  delta_X=0.052360, delta_Y=0.087266
# TILT MAGNITUDE: 0.101735
# TIME:  warping=19671.000000, total=22155
# COUNT: warping=1
```

The result is a space-separated list first containing the snapshot (first 6
columns) and current view (next 6 columns) image descriptor, followed by the
tilt correction parameters applied to the current view image (two columns:
delta_x and delta_y), the warping result (three columns: alpha, psi, d_min) and
finally a runtime measurement (three columns: warping count, warping time, total
time).

Following the space-separated results is a human-readable summary of the
results, prefixed by `#` to allow for filtering.

To evaluate the warping homing performance, compare the warping result (alpha,
psi and the derived home direction beta) to the true home direction, which can
be computed from the grid positions of the snapshot and current view descriptor
included in the output (the second and third columns following the filename).

Again, to get meaningful performance statistics, repeat and average this over
multiple image pairs and compare different tilt corrections by setting the
`tilt_solution` environment variable.

## TiltSearchDemo

The `TiltSearchDemo` application demonstrates the warping homing performance on
tilted panoramic images with searching for the best tilt correction parameters.
The application processes a single snapshot current view image pair with
optional orientation change. All other (warping and tilt correction related
parameters) are set via environment variables. The results are written to the
standard output:

```shell
# Usage: TiltSearchDemo ss_root ss cv_root cv [rotate-ss] [rotate-cv]
$ TiltSearchDemo data/lab/d 000000 data/lab/d 001000
000000 ... 001000 ... 0.02 0.14 0.147262 0.0490874 39136 225 16275 4166000 
# EST:   beta=3.043418, alpha=0.147262, psi=0.049087
# TRUE:  beta=3.141593, alpha=0.000000, psi=0.000000
# ERROR: beta=0.098175, alpha=0.147262, psi=0.049087
# TILT EST:   delta_X=0.020000, delta_Y=0.140000
# TILT TRUE:  delta_X=0.052360, delta_Y=0.087266
# TILT ERROR: 0.061760
# TIME:  warping=16275.000000, total=4166000
# COUNT: warping=225
```

The format of the results output is the same as for the `BaselineDemo`, but the
tilt correction parameters are now the found minimum and not given, ground truth
values.

The `TiltSearchDemo` application provides the additional option to produce a log
of the search process by setting the `logging` environment variable:

```shell
# Usage: TiltSearchDemo ss_root ss cv_root cv [rotate-ss] [rotate-cv]
$ search_strategy=1 logging=1 TiltSearchDemo data/lab/d 000000 data/lab/d 001000
000000 ... 001000 ... 0.035 0.105 0.0981748 0 39374 17 16592 320189 
> 0 0  0 0 41162 0.14 0.14  
> 0 0  0 0 41162 0.07 0.07  
> 0 0.07  0 0 40101 0.07 0.07  
> 0 0.07  0 0 40101 0.035 0.035  
> 0 0.105  0 0 39890 0.035 0.035  
> 0.035 0.105  0.0981748 0 39374 0.035 0.035  
> 0.035 0.105  0.0981748 0 39374 0.0175 0.0175  
# EST:   beta=3.043418, alpha=0.098175, psi=0.000000
# TRUE:  beta=3.141593, alpha=0.000000, psi=0.000000
# ERROR: beta=0.098175, alpha=0.098175, psi=0.000000
# TILT EST:   delta_X=0.035000, delta_Y=0.105000
# TILT TRUE:  delta_X=0.052360, delta_Y=0.087266
# TILT ERROR: 0.024760
# TIME:  warping=16592.000000, total=320189
# COUNT: warping=17
```

The log output follows the results and is prefixed by `>` to allow for
filtering. The actual format of the log depends on the chosen search strategy
and is intended for creating visualizations of the search process. The (per
line) format of the log's string representation is as follows:

#### Grid Searcher
Each line contains an N-dimensional coordinate in search space `xs` followed
by the objective function `value` at these coordinates.

```
# xs value
# here 2d + 3 warping values: delta_X delta_Y alpha psi d_min
-0.14 -0.14  6.18501 6.2341 44654
```

#### Pattern Searcher
Each line represents a single step of the search process containing the current
`argmin` (pattern center), the objective function value at this point `min` and
the `width` of the pattern as the distance from center per axis. The `argmin` is
an N-dimension coordinates in search space.

```
# argmin min width
# here 2d + 3 warping values: delta_X delta_Y alpha psi d_min width_X width_Y
0 0  0 0 41162 0.14 0.14
```

#### Nelder-Mead Searcher
Each line represents a single step of the search process containing the current
set of N+1 vertices, each followed by the objective function `value` at the
`vertex`. Each vertex is an N-dimension coordinate in search space.

```
# vertex value vertex value ...
# here 2d + 3 warping values: delta_X delta_Y alpha psi d_min ...
0 0.14  0.0981748 0.0490874 39425 -0.07 -0.07  6.2341 6.2341 42573 0.14 0  0.638136 0.245437 43507
```


### Plotting the Search Log
To run a tilt search and directly produce a plot visualizing the search process,
run the `searchplot` wrapper script:

```shell
# Usage: ./searchplot plot-file ss_root ss cv_root cv [rotate-ss] [rotate-cv]
$ search_strategy=1 searchplot pattern.svg data/lab/d 000000 data/lab/d 001000
```

Instead of printing the `TiltSearchDemo` results, this script filters the
output and forwards the logging to python scripts creating visualizations of the
search process. These scripts require `Python 3.5`, `numpy` and `matplotlib`.

## Environment Variables

To get a list of all environment variables read by one of the demo programs run
it in verbose mode (with the environment variable `verbose` set to `1`), e.g.
the `BaselineDemo`:

```shell
# Usage: BaselineDemo ss_root ss cv_root cv [rotate-ss] [rotate-cv] [no-tilt]
verbose=1 BaselineDemo data/lab/d 000000 data/lab/d 001000
# verbose = 1 (0) 
# dump_default = 1 (1) 
# nAlpha = 128 (128) 
# nPsi = 128 (128) 
# nScalePlanes = 9 (9) 
# ...
# tilt_solution = 0 (0) 
# interpolation = 0 (0) 
000000 ... 001000 ... 0.05236 0.087266 0.0981748 0 40246 1 19794 22316 
# EST:   beta=3.043418, alpha=0.098175, psi=0.000000
# TRUE:  beta=3.141593, alpha=0.000000, psi=0.000000
# ERROR: beta=0.098175, alpha=0.098175, psi=0.000000
# TILT EST:   delta_X=0.052360, delta_Y=0.087266
# TILT TRUE:  delta_X=0.052360, delta_Y=0.087266
# TILT MAGNITUDE: 0.101735
# TIME:  warping=19794.000000, total=22316
# COUNT: warping=1
```

Default values are shown in brackets. To list changed variables only, set
`dump_default=0`. Verbose output is prefixed by `#` so it can be ignored by
typical csv or table readers or filtered manually e.g. using `grep` or `sed`.

## C++17

The applications and the included tilt correction and search library make
extensive use of C++14 and C++17 features (return type deduction, class template
argument deduction, structured bindings, fold expression, ...), thus at least
`gcc 7` or `clang 5` are required to build the applications with `-std=c++17`.

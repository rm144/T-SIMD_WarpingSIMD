# Warping SIMD Tilt Correction+Search Library

If you are just interested in *using* the tilt correction search, have a look at
`TiltSearch.H` - this bundles all the requirements into a single, mostly runtime
configurable frontend. The template parameter `Model` selects the tilt parameter
model to use - either `RollPitch` or `AxisAngle`. The `Warping` parameter must
be a type of `WarpingBundle`, this especially selects the type of `SIMDImage`
to use and all other function (e.g. the tilt correction) are constructed to
match this. The following is a sketch of a simple warping tilt search using
default and hardcoded parameters - for more details and options see the
constructor documentation of the `tilt::Search` class.

```c++
// Construct a pattern searcher with
// - initial pattern centered at the origin
// - 0.14 radian wide along each axis
// - terminate if smaller than 0.02 radian on all axes
// - shrink by halving (x 0.5) the pattern
auto searcher = search::Pattern{
    std::tuple{0.0_rad, 0.0_rad}, 0.14, 0.02, 0.5
};
// Construct a box constraint function preventing the pattern from leaving the
// [-0.14;+0.14] radian range
auto limit = [X, Y](auto x, auto y) -> bool {
    return x >= -0.14_rad && x <= +0.14_rad && y >= -0.14_rad && y <= +0.14_rad;
};
// Construct a constrained tilt search of exact tilt correction parameterized by
// roll-pitch tilt given in radians (matching the searcher). Assumes warping is
// a previously constructed WarpingBundle.
auto tilt_search = tilt::Search{
    tilt::RollPitch<Radian>{}, tilt::EXACT, searcher, warping, limit
};
// Run the (constrained) tilt search on a snapshot-current view image pair;
// params is an instance of the WarpingBundle::Parameters structure.
auto[min, argmin] = tilt_search.run(ss, cv, params);
```

Full runtime configurability may be achieved by putting the construction (like
above) in a factory method handling different searchers via a `switch-case`
construct.

The following describes the components and utilities contained in this library
which are used to build and configure the tilt search frontend. Most of the
`Tilt*.H` components have dependencies on `SIMDImage.H` and
`SIMDImageFunctions.H`.

## Tilt Correction of Panoramic Images

Three solutions for tilt correction of panoramic images (this means
`SIMDImage`with `Panorama` addon) are implemented: `EXACT`, `APPROXIMATE`,
`VERTICAL`. These are implemented in `TiltCorrection.H` as overloads of the
`tilt::correction` function template. The image transformation is done using
the `remap` and `shift` functions from `SIMDImageFunctions.H`. At compile time
selection of the tilt solution is most easily done via the first template
argument of the `tilt::correct` hub template, this selects the concrete
implementation via tag dispatching. A runtime selection is possible via the
`tilt::make_correction` factory function returning an instance of the
`tilt::Correction` interface type wrapping the implementation in a runtime
polymorphic `std::function`. While this is more flexible by allowing
construction via a runtime provided selector number, it may introduce a *slight*
runtime overhead. The runtime selector numbers are:

    EXACT=0, APPROXIMATE=1, VERTICAL=2

The tilt correction functions expect the image tilt parameterized as axis-angle
representation (rather angle-to-axis - angle), while most of *the outside* (e.g.
the databases) represents tilt as roll-pitch. Both representations are
implemented as so called *tilt models* (`RollPitch` and `AxisAngle`) in
`TiltModel.H` with implicit conversion operators between them. Thus, it is
perfectly fine to pass `RollPitch` parameterized tilt to the correction function
expecting `AxisAngle` - the appropriate conversion is automatically inserted by
the compiler.

## Generic Search Strategies

Headers starting with `Search*.H` contain structures, functions and utilities
for optimization (~ minimization) of generic objective functions. The frontend
of the generic search is the `search::minimize` function template and related
functions. Once a search method provides a minimization function, maximization
and constrained optimization are derived automatically.

### Objective Function and Type Requirements

A suitable objective function is any function mapping from a set of parameter
values to a single scalar-like objective value. Typically, this would be a
function object like a lambda expression which is allowed to capture state (even
by reference) as long as this does not go out of scope before the search
terminates. The objective function might even be a template in some
circumstances, if it is possible to unambiguously deduce the types from the
searcher (and optionally constraint function) passed to the optimizer.

```c++
// Construct a quadratic test function with minimum of c0 at (x0,y0)
auto f =[x0, y0, c0](const double x, const double y) {
    return (x - x0) * (x - x0) + (y - y0) * (y - y0) + c0;
};
// Run an unconstrained search with some previously configured searcher to find
// the minimum
auto[min, argmin] = search::minimize(f, searcher);
```

The parameter values may be of various, even differing, types but must behave
roughly like a floating point type: this means, it must be default
constructible, copyable and arithmetic, (compound) assignment and comparison
operators must be available. For the objective value it is sufficient to have a
default constructible, copyable and comparable type for *most* of the search
strategies. Only the Nelder-Mead searcher, when using the *sample standard
deviation of the function values* as the termination criterion, requires
arithmetic and arithmetic assignment operators for the objective value type.

For constrained optimization, a reasonable value for *infinity* is required. 
However, instead of specializing the standard C++ class template
`std::numeric_limits<Type>::infinity()`, it may be sufficient to provide the
infinity value by default initialization just like the `WarpingBundle::Result`
type.

If the logging functionality from `SearchLogger.H` is used, it is in addition
necessary to provide overloads of output stream (`std::ostream`) insertion
operators for parameter and objective value types to build the log's string
representation - in case of builtin types like `float`/`double` this is already
the case.

### Constrained Optimization

The optimization methods accept an optional third parameter for specifying a
constraint function. This function must take the same parameters as the
objective function but returns `true` or `false` on valid or invalid inputs
respectively. The constraint function is evaluated first - before evaluating the
objective - and the objective function is not called at all, if it yields
`false`.

### Search Methods

Three direct search methods (no gradients required, only doing comparisons) are
currently implemented: the exhaustive grid search and the heuristic pattern and
Nelder-Mead (aka downhill simplex) search. All search methods are designed to
work with `search::minimize` or its constrained variant. The searchers do not do
the full search on their own, but merely operate on the search space and
encapsulate state via the `init`, `step`, `done` or `for_each` methods.

When swapping the searcher, no modifications to the objective or call site are
necessary - except for passing another searcher. However, it is not directly
possible to select or swap the searcher at runtime as `search::minimize` is a
pure static, compile time polymorphic interface. Runtime selection can be done
via an extra indirection in form of a runtime polymorphic wrapper using virtual
calls (see `TiltSearch.H` for an example).

Note: The search methods differ vastly in their construction, initial and
termination parameters! See the constructor documentation of `NelderMead` or
`Pattern` for details.

### Search Process Logging

A wrapper to add logging of `std::string` representations of the search process
into a `std::list` is provided in `SearchLogger.H`. This list can then be
printed or saved to file. In case of the grid searcher, this overloads the
`search::minimize` function, in case of the heuristic searchers, it hooks the
`init` and `step` methods.

```c++
// Set up an empty list to write the log into
std::list<std::string> log;
// Wrap the already configured searcher with the logger writing to the list
// passed by reference and run the search to find the minimum
auto[min, argmin] = search::minimize(f, search::Logger(searcher, log));
// Print the log one state per line
for (const auto& state : log) {
    std::cout << state << std::endl;
}
```

As logging requires memory allocations and formatted string stream insertions,
it is not recommended doing logging and performance measurements at the same
time. The (per line) format of the log's string representation is as follows:

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

## Utilities
To avoid angular unit confusions, `Angle.H` provides a strong type angle type
wrapping basic types. The tilt models and therefor the tilt corrections make
extensive use of this. Via implicit conversion operators it is still possible to
e.g. pass degree to a function excepting radian (the compiler will generate the
conversion automatically), but construction and value extraction (e.g. back to a
plain `double`) must be made explicit.

The pattern and Nelder-Mead searcher need to do arithmetic and comparison
operations on `std::tuple` of different types. These operations - e.g.
elementwise addition of two tuples - are implemented in `TupleUtility.H` in
their own namespace. This especially includes `std::ostream` insertion of
`std::tuple` necessary to build the search log string representations.

## C++17
This library makes extensive use of C++14 and C++17 features (return type
deduction, class template argument deduction, structured bindings, fold
expression, ...), thus at least `gcc 7` or `clang 5` are required to build
applications including from this library with `-std=c++17` set.

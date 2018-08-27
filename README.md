# ImplicitFunctions.jl

Implements a simple Iterator to track continuous branches of implicit functions
by combining ForwardDiff.jl and NLsolve.jl.

Input functions should be of the form <img src="/tex/e363ebb8835aa265961720caf0517f78.svg?invert_in_darkmode&sanitize=true" align=middle width=77.43130889999999pt height=24.65753399999998pt/>. Given a starting
point <img src="/tex/a6ee443b8b9ffeb5de0c462c42d96036.svg?invert_in_darkmode&sanitize=true" align=middle width=53.61303089999999pt height=24.65753399999998pt/> the Iterator then produces a series of points on
<img src="/tex/f98a972dbd244305c46ceda2aac5c780.svg?invert_in_darkmode&sanitize=true" align=middle width=31.558228349999993pt height=24.65753399999998pt/> by first making a predictive step along the tangent of the
implicit curve and then correcting it using Newton's method.

The current implementation is very naive and likely produces some unexpected
behavior; especially at bifurcation points. Furthermore, it requires functions
implemented in pure julia, in order to calculate Jacobians using ForwardDiff.jl.
Currently, the algorithm is fixed-step-only (in <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/>) which means only branches
monotonical in <img src="/tex/2ec6e630f199f589a2402fdf3e0289d5.svg?invert_in_darkmode&sanitize=true" align=middle width=8.270567249999992pt height=14.15524440000002pt/> can be mapped.

## Usage

Specify the function <img src="/tex/3cf4fbd05970446973fc3d9fa3fe3c41.svg?invert_in_darkmode&sanitize=true" align=middle width=8.430376349999989pt height=14.15524440000002pt/>, the starting point <img src="/tex/a6ee443b8b9ffeb5de0c462c42d96036.svg?invert_in_darkmode&sanitize=true" align=middle width=53.61303089999999pt height=24.65753399999998pt/> and the
step-width <img src="/tex/4803b03804422d02acec2246c34613f8.svg?invert_in_darkmode&sanitize=true" align=middle width=16.19863904999999pt height=22.831056599999986pt/>:

```julia
g(u::AbstractVector, p:::Real) = ...
u0 = zeros(...)
p0 = 1.0
dp = 0.1
iter = ImplicitIterator(g, u0, p0, dp)

# get the first 5 points of the implicit function
collect(Iterators.take(iter, 5))
```

To get debugging information use:

```julia
iter = ImplicitIterator(g, u0, p0, dp, debug = true)
```

## Versions

Please note, because the Iterator interface has changed between julia version
0.6 and 0.7 please use the [julia-0.6](https://github.com/gwater/ImplicitFunctions.jl/tree/julia-0.6) branch for older projects.

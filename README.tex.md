# ImplicitFunctions.jl

Implements a simple Iterator to track continuous branches of implicit functions
by combining ForwardDiff.jl and NLsolve.jl.

Input functions should be of the form $g:\mathbb{R}^n \times \mathbb{R}
\rightarrow \mathbb{R}^n$. Given a starting point $(\mathbf{u}_0,p_0)$ the
Iterator then produces a series of points on $\mathbf{u}(p)$ such that
everywhere $g(\mathbf{u},p)=0$ by first making a predictive step along the
tangent of the implicit curve and then correcting it using Newton's method.

The current implementation is very naive and likely produces some unexpected
behavior; especially at bifurcation points. Furthermore, it requires functions
implemented in pure julia in order to calculate Jacobians using ForwardDiff.jl.
Currently, the algorithm is fixed-step-only (in $p$) which means only branches
with strictly increasing / decreasing $p$ can be mapped.

## Usage

Specify the function $g$, the starting point $(\mathbf{u}_0,p_0)$ and the
step-width $\delta p$:

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

## See also

Also check out [ImplicitEquations.jl](https://github.com/jverzani/ImplicitEquations.jl) which uses different numerical approach to graphically display implicitly defined regions.

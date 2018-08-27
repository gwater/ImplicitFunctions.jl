module ImplicitFunctions

using ForwardDiff
using Tensors
using NLsolve

import Base: start, next, done, iteratorsize, eltype

export ImplicitIterator
export start, next, done, iteratorsize, eltype

jacobian(f, x::AbstractVector) =
    Tensors.gradient(x -> Vec{length(x)}(f(x)), Vec{length(x)}(x))
jacobian(f, x::Real) = ForwardDiff.derivative(f, x)

# https://math.stackexchange.com/questions/1415898/pseudo-arclength-continuation-scheme
function arc_approx(f, x0, p0::P, dp::P) where {P <: Real}
    # tangent vector
    dxdp = -jacobian(x -> f(x, p0), x0) \ jacobian(p -> f(x0, p), p0)
    return x0 .+ dxdp * dp, p0 + dp
end

struct ConvergenceException{T} <: Exception
    res::T
end

function _newton(f, x0, ::Type{Val{DEBUG}}) where DEBUG
    res = nlsolve(f, x0, method = :newton, inplace = false, autodiff = :forward)
    DEBUG && info(res)
    converged(res) && return res.zero
    throw(ConvergenceException(res))
end

struct ImplicitIterator{X, P <: Real, DEBUG, F}
    f::F
    x0::X
    p0::P
    dp::P
end
ImplicitIterator(f::F, x0::X, p0::P, dp::P = 0.1; debug::Bool = false) where
    {X, P, F} = ImplicitIterator{X, P, debug, F}(f, x0, p0, dp)

iteratorsize(::Type{I}) where {I <: ImplicitIterator} = Base.SizeUnknown()
eltype(::Type{I}) where {X, P, I <: ImplicitIterator{X, P}} = Tuple{X, P}

function start(iter::ImplicitIterator)
    return iter.x0, iter.p0
end

function next(iter::ImplicitIterator{X, P, DEBUG}, old_state) where {X, P, DEBUG}
    x0, p0 = old_state
    x1_approx, p1 = arc_approx(iter.f, x0, p0, iter.dp)
    x1 = _newton(x -> iter.f(x, p1), x1_approx, Val{DEBUG})
    new_state = x1, p1
    return new_state, new_state
end

function done(iter::ImplicitIterator{X, P, DEBUG}, state) where {X, P, DEBUG}
    # stop if we cannot invert the Jacobian df/dx
    x0, p0 = state
    J = jacobian(x -> iter.f(x, p0), x0)
    DEBUG && info(J)
    return iszero(det(J)) # break at non-invertible Jacobian
end

end # module

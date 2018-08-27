module ImplicitFunctions

using ForwardDiff
using Tensors
using NLsolve

import Base: iterate, eltype, IteratorSize

export ImplicitIterator
export iterate, eltype, IteratorSize

jacobian(f, x::AbstractVector) =
    Tensors.gradient(x -> Vec{length(x)}(f(x)), Vec{length(x)}(x))
jacobian(f, x::Real) = ForwardDiff.derivative(f, x)

# https://math.stackexchange.com/questions/1415898/pseudo-arclength-continuation-scheme
function arc_approx(f, x0, p0::P, dp::P, ::Type{Val{DEBUG}} = Val{false}) where {P <: Real, DEBUG}
    # tangent vector
    Jx = jacobian(x -> f(x, p0), x0)
    DEBUG && info(Jx)
    # iszero(det(Jx)) % throw some exception
    dxdp = -Jx \ jacobian(p -> f(x0, p), p0)
    return x0 .+ dxdp * dp, p0 + dp
end

struct ConvergenceException{T} <: Exception
    res::T
end

function _newton(f, x0, ::Type{Val{DEBUG}} = Val{false}) where DEBUG
    res = nlsolve(f, x0, method = :newton, inplace = false, autodiff = :forward)
    DEBUG && info(res)
    converged(res) && return res.zero
    throw(ConvergenceException(res))
end

function step(f, x0, p0, dp, ::Type{Val{DEBUG}} = Val{false}) where DEBUG
    x1_approx, p1 = arc_approx(f, x0, p0, dp, Val{DEBUG})
    x1 = _newton(x -> f(x, p1), x1_approx, Val{DEBUG})
    return x1, p1
end

# Iterator interface
# https://docs.julialang.org/en/v1.0.0/manual/interfaces/#man-interface-iteration-1

struct ImplicitIterator{X, P <: Real, DEBUG, F}
    f::F
    x0::X
    p0::P
    dp::P
end
ImplicitIterator(f::F, x0::X, p0::P, dp::P = 0.1; debug::Bool = false) where
    {X, P, F} = ImplicitIterator{X, P, debug, F}(f, x0, p0, dp)

IteratorSize(::Type{I}) where {I <: ImplicitIterator} = Base.SizeUnknown()

eltype(::Type{I}) where {X, P, I <: ImplicitIterator{X, P}} = Tuple{X, P}

function iterate(iter::ImplicitIterator)
    initial_state = (iter.x0, iter.p0)
    return initial_state, initial_state
end

function iterate(iter::ImplicitIterator{X, P, DEBUG}, old_state) where {X, P, DEBUG}
    x0, p0 = old_state
    new_state = step(iter.f, x0, p0, iter.dp, Val{DEBUG})
    return new_state, new_state
end

end # module

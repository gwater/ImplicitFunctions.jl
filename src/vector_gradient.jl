module VectorGradient

using StaticArrays
import ForwardDiff: Dual, partials

export gradient

@generated function _partial_dual(v::T, ::Type{Val{N}}, ::Type{Val{I}}) where {T, N, I}
    zs = zeros(T, N)
    zs[I] = one(T)
    return quote
        $(Expr(:meta, :inline))
        return Dual(v, $(zs...))
    end
end

@inline function _wrap_dual(v::SVector{N, T}) where {N, T}
    DT = Dual{Nothing, T, N}
    @inbounds v_dual = SVector{N, DT}(NTuple{N, DT}(_partial_dual(v[i], Val{N}, Val{i}) for i in 1:N))
    return v_dual
end

@generated function _extract_jacobian(v::SVector{M}, ::SVector{N}) where {M, N}
    pp = [:(partials(v[$i])) for i in 1:M]
    return quote
        $(Expr(:meta, :inline))
        return SMatrix{M, N}(vcat($(pp...)))
    end
end

function gradient(f::F, v::V) where {F, N, V <: SVector{N}}
    v_dual = _wrap_dual(v)
    res = f(v_dual)
    return _extract_jacobian(res, v)
end

end #module

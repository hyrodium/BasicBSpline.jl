# B-spline manifold

"""
B-spline manifold for lower polynomial degree
TODO: make the field `bsplinespaces` to be conposite type, not abstract type, for performance
"""
struct FastBSplineManifold <: AbstractBSplineManifold
    # TODO: avoid abstract type FastBSplineSpace
    bsplinespaces::Array{T,1} where T <: FastBSplineSpace
    controlpoints::Array{Float64}
    function FastBSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::AbstractArray{<:Real})
        Ps = FastBSplineSpace.(Ps)
        if collect(size(a)[1:end-1]) ≠ dim.(Ps)
            error("dimension does not match")
        else
            P = convert(Array{FastBSplineSpace,1}, Ps)
            a = convert(Array{Float64}, a)
            new(P, a)
        end
    end
end

function FastBSplineManifold(Ps::AbstractArray{<:AbstractBSplineSpace,1}, a::Array{Array{T,1}} where T<:Real)
    d̂ = length(a[1])
    A = reshape(transpose(hcat(reshape(a,prod(size(a)))...)), size(a)..., d̂)
    return FastBSplineManifold(Ps, A)
end

"""
convert AbstractBSplineManifold to FastBSplineManifold
"""
function FastBSplineManifold(M::AbstractBSplineManifold)
    FastBSplineManifold(BSplineSpace.(M.bsplinespaces), M.controlpoints)
end

@doc raw"""
Calculate the mapping of B-spline manifold for given parameter.
```math
\bm{p}(t^1,\dots,t^d)
=\sum_{i^1,\dots,i^d}B_{i^1,\dots,i^d}(t^1,\dots,t^d) \bm{a}_{i^1,\dots,i^d}
```
"""
function mapping(M::FastBSplineManifold, t::Array{<:Real,1})
    # TODO: faster
    return mapping(BSplineManifold(M),t)
end

@doc raw"""
Calculate the dimension of B-spline manifold.
"""
function dim(M::FastBSplineManifold)
    length(M.bsplinespaces)
end

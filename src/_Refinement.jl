# Refinement

# TODO: general dimension
# TODO: Update docstrings

@doc raw"""
Refinement of B-spline manifold with given B-spline spaces.
"""
refinement

function refinement(M::BSplineManifold{1}, Ps′::Tuple{<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,) = dim.(Ps)
    (n1′,) = dim.(Ps′)
    (A1,) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * a[I₁] for I₁ in 1:n1) for J₁ in 1:n1′]
    return BSplineManifold(a′, Ps′)
end
function refinement(M::BSplineManifold{2}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,n2) = dim.(Ps)
    (n1′,n2′) = dim.(Ps′)
    (A1,A2) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′]
    return BSplineManifold(a′, Ps′)
end
function refinement(M::BSplineManifold{3}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,n2,n3) = dim.(Ps)
    (n1′,n2′,n3′) = dim.(Ps′)
    (A1,A2,A3) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * a[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′]
    return BSplineManifold(a′, Ps′)
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knotvector.
"""
function refinement(M::BSplineManifold; p₊::Union{Nothing,NTuple{Dim,Int}}=nothing, k₊::Union{Nothing,NTuple{Dim,KnotVector{T}}}=nothing) where {Dim, T}
    Ps = collect(bsplinespaces(M))
    d = length(Ps)
    if isnothing(p₊)
        p₊ = zeros(Int, d)
    elseif length(Ps) ≠ length(p₊)
        throw(DimensionMismatch())
    end
    if isnothing(k₊)
        k₊ = zeros(KnotVector, d)
    elseif length(Ps) ≠ length(k₊)
        throw(DimensionMismatch())
    end

    Ps′ = [expandspace(Ps[i], p₊=p₊[i], k₊=k₊[i]) for i in 1:Dim]

    return refinement(M, tuple(Ps′...))
end

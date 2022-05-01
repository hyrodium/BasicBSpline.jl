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

function refinement(M::RationalBSplineManifold{1}, Ps′::Tuple{<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    (n1,) = dim.(Ps)
    (n1′,) = dim.(Ps′)
    (A1,) = changebasis.(Ps, Ps′)

    w′ = [sum(A1[I₁,J₁] * w[I₁] for I₁ in 1:n1) for J₁ in 1:n1′]
    a′ = [sum(A1[I₁,J₁] * a[I₁] * w[I₁] for I₁ in 1:n1) for J₁ in 1:n1′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end
function refinement(M::RationalBSplineManifold{2}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    (n1,n2) = dim.(Ps)
    (n1′,n2′) = dim.(Ps′)
    (A1,A2) = changebasis.(Ps, Ps′)

    w′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * w[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′]
    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂] * w[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end
function refinement(M::RationalBSplineManifold{3}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    w = weights(M)
    (n1,n2,n3) = dim.(Ps)
    (n1′,n2′,n3′) = dim.(Ps′)
    (A1,A2,A3) = changebasis.(Ps, Ps′)

    w′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * w[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′]
    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * a[I₁,I₂,I₃] * w[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knotvector.
"""
function refinement(M::AbstractManifold{Dim}; p₊::Union{Nothing,NTuple{Dim,Int}}=nothing, k₊::Union{Nothing,NTuple{Dim,KnotVector{T} where T}}=nothing) where Dim
    Ps = bsplinespaces(M)
    if isnothing(p₊) & isnothing(k₊)
        Ps′ = Ps
    elseif isnothing(p₊) & !isnothing(k₊)
        Ps′ = [expandspace(Ps[i], k₊=k₊[i]) for i in 1:Dim]
    elseif !isnothing(p₊) & isnothing(k₊)
        Ps′ = [expandspace(Ps[i], p₊=p₊[i]) for i in 1:Dim]
    else
        Ps′ = [expandspace(Ps[i], p₊=p₊[i], k₊=k₊[i]) for i in 1:Dim]
    end

    return refinement(M, tuple(Ps′...))
end

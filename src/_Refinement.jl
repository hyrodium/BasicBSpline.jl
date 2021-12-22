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
    d = size(a,2)  # 2 == Dim+1
    (n1,) = dim.(Ps)
    (n1′,) = dim.(Ps′)
    (A1,) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * a[I₁,j] for I₁ in 1:n1) for J₁ in 1:n1′, j in 1:d]
    return BSplineManifold(Ps′, a′)
end
function refinement(M::BSplineManifold{2}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    d = size(a,3)  # 3 == Dim+1
    (n1,n2) = dim.(Ps)
    (n1′,n2′) = dim.(Ps′)
    (A1,A2) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂,j] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′, j in 1:d]
    return BSplineManifold(Ps′, a′)
end
function refinement(M::BSplineManifold{3}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    d = size(a,4)  # 4 == Dim+1
    (n1,n2,n3) = dim.(Ps)
    (n1′,n2′,n3′) = dim.(Ps′)
    (A1,A2,A3) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂]* A3[I₃,J₃] * a[I₁,I₂,I₃,j] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′, j in 1:d]
    return BSplineManifold(Ps′, a′)
end

function refinement(M::CustomBSplineManifold{1}, Ps′::Tuple{<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,) = dim.(Ps)
    (n1′,) = dim.(Ps′)
    (A1,) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * a[I₁] for I₁ in 1:n1) for J₁ in 1:n1′]
    return CustomBSplineManifold(Ps′, a′)
end
function refinement(M::CustomBSplineManifold{2}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,n2) = dim.(Ps)
    (n1′,n2′) = dim.(Ps′)
    (A1,A2) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′]
    return CustomBSplineManifold(Ps′, a′)
end
function refinement(M::CustomBSplineManifold{3}, Ps′::Tuple{<:AbstractBSplineSpace,<:AbstractBSplineSpace,<:AbstractBSplineSpace})
    Ps = bsplinespaces(M)
    a = controlpoints(M)
    (n1,n2,n3) = dim.(Ps)
    (n1′,n2′,n3′) = dim.(Ps′)
    (A1,A2,A3) = changebasis.(Ps, Ps′)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂]* A3[I₃,J₃] * a[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′]
    return CustomBSplineManifold(Ps′, a′)
end


@doc raw"""
Refinement of B-spline manifold with additional degree and knots.
"""
function refinement(M::BSplineManifold{Dim}; p₊::Union{Nothing,NTuple{Dim,Int}}=nothing, k₊::Union{Nothing,NTuple{Dim,Knots{T}}}=nothing) where {Dim, T}
    Ps = bsplinespaces(M)
    if isnothing(p₊)
        p₊ = tuple(zeros(Int, Dim)...)
    elseif length(Ps) ≠ length(p₊)
        throw(DimensionMismatch())
    end
    if isnothing(k₊)
        k₊ = tuple(zeros(Knots, Dim)...)
    elseif length(Ps) ≠ length(k₊)
        throw(DimensionMismatch())
    end

    Ps′ = [(P = Ps[i];
    p = degree(P);
    k = knots(P);
    k_unique = unique(k[1+p:end-p]);
    BSplineSpace{p + p₊[i]}(k + p₊[i] * k_unique + k₊[i])) for i in 1:Dim]

    return refinement(M, tuple(Ps′...))
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knots.
"""
function refinement(M::CustomBSplineManifold; p₊::Union{Nothing,NTuple{Dim,Int}}=nothing, k₊::Union{Nothing,NTuple{Dim,Knots{T}}}=nothing) where {Dim, T}
    Ps = collect(bsplinespaces(M))
    𝒂 = controlpoints(M)
    d = length(Ps)
    n = dim.(Ps)
    if isnothing(p₊)
        p₊ = zeros(Int, d)
    elseif length(Ps) ≠ length(p₊)
        throw(DimensionMismatch())
    end
    if isnothing(k₊)
        k₊ = zeros(Knots, d)
    elseif length(Ps) ≠ length(k₊)
        throw(DimensionMismatch())
    end

    Ps′ = [(P = Ps[i];
    p = degree(P);
    k = knots(P);
    k_unique = unique(k[1+p:end-p]);
    BSplineSpace{p + p₊[i]}(k + p₊[i] * k_unique + k₊[i])) for i in 1:Dim]

    return refinement(M, tuple(Ps′...))
end

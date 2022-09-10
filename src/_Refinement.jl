# Refinement

# TODO: general dimension
# TODO: Update docstrings

@doc raw"""
Refinement of B-spline manifold with given B-spline spaces.
"""
refinement

function refinement(M::BSplineManifold{1}, Ps′::NTuple{1, AbstractBSplineSpace})
    P1, = bsplinespaces(M)
    P1′, = Ps′
    n1 = dim(P1)
    n1′ = dim(P1′)
    A1 = changebasis(P1, P1′)
    a = controlpoints(M)

    a′ = [sum(A1[I₁,J₁] * a[I₁] for I₁ in 1:n1) for J₁ in 1:n1′]
    return BSplineManifold(a′, Ps′)
end
function refinement(M::BSplineManifold{2}, Ps′::NTuple{2, AbstractBSplineSpace})
    P1, P2 = bsplinespaces(M)
    P1′, P2′ = Ps′
    n1 = dim(P1)
    n2 = dim(P2)
    n1′ = dim(P1′)
    n2′ = dim(P2′)
    A1 = changebasis(P1, P1′)
    A2 = changebasis(P2, P2′)
    a = controlpoints(M)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′]
    return BSplineManifold(a′, Ps′)
end
function refinement(M::BSplineManifold{3}, Ps′::NTuple{3, AbstractBSplineSpace})
    P1, P2, P3 = bsplinespaces(M)
    P1′, P2′, P3′ = Ps′
    n1 = dim(P1)
    n2 = dim(P2)
    n3 = dim(P3)
    n1′ = dim(P1′)
    n2′ = dim(P2′)
    n3′ = dim(P3′)
    A1 = changebasis(P1, P1′)
    A2 = changebasis(P2, P2′)
    A3 = changebasis(P3, P3′)
    a = controlpoints(M)

    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * a[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′]
    return BSplineManifold(a′, Ps′)
end

function refinement(M::RationalBSplineManifold{1}, Ps′::NTuple{1, AbstractBSplineSpace})
    P1, = bsplinespaces(M)
    P1′, = Ps′
    n1 = dim(P1)
    n1′ = dim(P1′)
    A1 = changebasis(P1, P1′)
    a = controlpoints(M)
    w = weights(M)

    w′ = [sum(A1[I₁,J₁] * w[I₁] for I₁ in 1:n1) for J₁ in 1:n1′]
    a′ = [sum(A1[I₁,J₁] * a[I₁] * w[I₁] for I₁ in 1:n1) for J₁ in 1:n1′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end
function refinement(M::RationalBSplineManifold{2}, Ps′::NTuple{2, AbstractBSplineSpace})
    P1, P2 = bsplinespaces(M)
    P1′, P2′ = Ps′
    n1 = dim(P1)
    n2 = dim(P2)
    n1′ = dim(P1′)
    n2′ = dim(P2′)
    A1 = changebasis(P1, P1′)
    A2 = changebasis(P2, P2′)
    a = controlpoints(M)
    w = weights(M)

    w′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * w[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′]
    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * a[I₁,I₂] * w[I₁,I₂] for I₁ in 1:n1, I₂ in 1:n2) for J₁ in 1:n1′, J₂ in 1:n2′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end
function refinement(M::RationalBSplineManifold{3}, Ps′::NTuple{3, AbstractBSplineSpace})
    P1, P2, P3 = bsplinespaces(M)
    P1′, P2′, P3′ = Ps′
    n1 = dim(P1)
    n2 = dim(P2)
    n3 = dim(P3)
    n1′ = dim(P1′)
    n2′ = dim(P2′)
    n3′ = dim(P3′)
    A1 = changebasis(P1, P1′)
    A2 = changebasis(P2, P2′)
    A3 = changebasis(P3, P3′)
    a = controlpoints(M)
    w = weights(M)

    w′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * w[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′]
    a′ = [sum(A1[I₁,J₁] * A2[I₂,J₂] * A3[I₃,J₃] * a[I₁,I₂,I₃] * w[I₁,I₂,I₃] for I₁ in 1:n1, I₂ in 1:n2, I₃ in 1:n3) for J₁ in 1:n1′, J₂ in 1:n2′, J₃ in 1:n3′] ./ w′
    return RationalBSplineManifold(a′, w′, Ps′)
end

function refinement(M::AbstractManifold{Dim}, Ps′::Vararg{AbstractBSplineSpace, Dim}) where Dim
    return refinement(M, Ps′)
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knotvector.
"""
@generated function refinement(M::AbstractManifold{Dim},
                               p₊::NTuple{Dim, Val},
                               k₊::NTuple{Dim, AbstractKnotVector}=ntuple(i->EmptyKnotVector(), Val(Dim))) where Dim
    Ps = [Symbol(:P,i) for i in 1:Dim]
    Ps′ = [Symbol(:P,i,"′") for i in 1:Dim]
    ks = [Symbol(:k,i,:₊) for i in 1:Dim]
    ps = [Symbol(:p,i,:₊) for i in 1:Dim]
    exP = Expr(:tuple, Ps...)
    exP′ = Expr(:tuple, Ps′...)
    exk = Expr(:tuple, ks...)
    exp = Expr(:tuple, ps...)
    exs = [:($(Symbol(:P,i,"′")) = expandspace($(Symbol(:P,i)), $(Symbol(:p,i,:₊)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exp = p₊),
        :($exk = k₊),
        exs...,
        :(return refinement(M, $(exP′)))
    )
end

@generated function refinement(M::AbstractManifold{Dim},
                               k₊::NTuple{Dim, AbstractKnotVector}=ntuple(i->EmptyKnotVector(), Val(Dim))) where Dim
    Ps = [Symbol(:P,i) for i in 1:Dim]
    Ps′ = [Symbol(:P,i,"′") for i in 1:Dim]
    ks = [Symbol(:k,i,:₊) for i in 1:Dim]
    exP = Expr(:tuple, Ps...)
    exP′ = Expr(:tuple, Ps′...)
    exk = Expr(:tuple, ks...)
    exs = [:($(Symbol(:P,i,"′")) = expandspace($(Symbol(:P,i)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exk = k₊),
        exs...,
        :(return refinement(M, $(exP′)))
    )
end

# resolve ambiguities
refinement(M::BasicBSpline.AbstractManifold{0}, ::Tuple{}) = M

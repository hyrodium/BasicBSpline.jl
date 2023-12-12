# Refinement

# TODO: general dimension
# TODO: Update docstrings

@doc raw"""
Refinement of B-spline manifold with given B-spline spaces.
"""
refinement

function _i_ranges_R(A, P′)
    n, n′ = size(A)
    i_ranges = fill(1:0, n′)
    i = 1
    for j in 1:n′
        isdegenerate_R(P′, j) && continue
        iszero(view(A,:,j)) && continue
        while true
            if Base.isstored(A,i,j)
                i_ranges[j] = i:n
                break
            end
            i = i + 1
        end
    end
    i = n
    for j in reverse(1:n′)
        isdegenerate_R(P′, j) && continue
        iszero(view(A,:,j)) && continue
        while true
            if Base.isstored(A,i,j)
                i_ranges[j] = first(i_ranges[j]):i
                break
            end
            i = i - 1
        end
    end
    return i_ranges
end

function _i_ranges_I(A, P′)
    n, n′ = size(A)
    i_ranges = fill(1:0, n′)
    i = 1
    for j in 1:n′
        isdegenerate_I(P′, j) && continue
        while true
            if Base.isstored(A,i,j)
                i_ranges[j] = i:n
                break
            end
            i = i + 1
        end
    end
    i = n
    for j in reverse(1:n′)
        isdegenerate_I(P′, j) && continue
        while true
            if Base.isstored(A,i,j)
                i_ranges[j] = first(i_ranges[j]):i
                break
            end
            i = i - 1
        end
    end
    return i_ranges
end

const _i_ranges = _i_ranges_R

function _ref_ctrl_elm(a::Array{T, Dim}, A::NTuple{Dim, SparseMatrixCSC}, R::NTuple{Dim, Vector{<:UnitRange}}, J::CartesianIndex{Dim}) where {T, Dim}
    ci = CartesianIndices(getindex.(R, J.I))
    if isempty(ci)
        # Should be type-stable
        S = Base.promote_op(*, eltype.(A)..., T)
        return zero(S)
    else
        return sum(prod(getindex.(A, I.I, J.I)) * a[I] for I in ci)
    end
end

function refinement_R(M::BSplineManifold{Dim}, P′::NTuple{Dim, BSplineSpace}) where Dim
    A = changebasis_R.(bsplinespaces(M), P′)
    R = _i_ranges_R.(A, P′)
    a = controlpoints(M)
    a′ = [_ref_ctrl_elm(a,A,R,J) for J in CartesianIndices(UnitRange.(1, dim.(P′)))]
    return BSplineManifold(a′, P′)
end

function refinement_I(M::BSplineManifold{Dim}, P′::NTuple{Dim, BSplineSpace}) where Dim
    A = changebasis_I.(bsplinespaces(M), P′)
    R = _i_ranges_I.(A, P′)
    a = controlpoints(M)
    a′ = [_ref_ctrl_elm(a,A,R,J) for J in CartesianIndices(UnitRange.(1, dim.(P′)))]
    return BSplineManifold(a′, P′)
end

function refinement(M::BSplineManifold{Dim}, P′::NTuple{Dim, BSplineSpace}) where Dim
    A = changebasis.(bsplinespaces(M), P′)
    R = _i_ranges.(A, P′)
    a = controlpoints(M)
    a′ = [_ref_ctrl_elm(a,A,R,J) for J in CartesianIndices(UnitRange.(1, dim.(P′)))]
    return BSplineManifold(a′, P′)
end

function refinement(M::RationalBSplineManifold{1}, Ps′::NTuple{1, BSplineSpace})
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
function refinement(M::RationalBSplineManifold{2}, Ps′::NTuple{2, BSplineSpace})
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
function refinement(M::RationalBSplineManifold{3}, Ps′::NTuple{3, BSplineSpace})
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

function refinement(M::AbstractManifold{Dim}, Ps′::Vararg{BSplineSpace, Dim}) where Dim
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
refinement(M::AbstractManifold{0}, ::Tuple{}) = M
refinement(M::BSplineManifold{0, Deg, C, S} where {Deg, C, S<:Tuple{}}, ::Tuple{}) = M

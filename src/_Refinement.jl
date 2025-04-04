# Refinement

# TODO: Update docstrings

@doc raw"""
Refinement of B-spline manifold with given B-spline spaces.
"""
refinement

function _i_ranges_R(A::AbstractSparseMatrix, P′::BSplineSpace)
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

function _i_ranges_I(A::AbstractSparseMatrix, P′::BSplineSpace)
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

function __control_points(value::T, index::Integer, dims::NTuple{Dim, Integer}) where {T, Dim}
    a = Array{T,Dim}(undef, dims)
    i = prod(index)
    for i in 1:(i-1)
        @inbounds a[i] = zero(value)
    end
    @inbounds a[i] = value
    return a
end

function _isempty(R::NTuple{1, Vector{UnitRange{Int}}}, J::CartesianIndex{1})
    return isempty(R[1][J[1]])
end
function _isempty(R::NTuple{2, Vector{UnitRange{Int}}}, J::CartesianIndex{2})
    return isempty(R[1][J[1]]) || isempty(R[2][J[2]])
end
function _isempty(R::NTuple{3, Vector{UnitRange{Int}}}, J::CartesianIndex{3})
    return isempty(R[1][J[1]]) || isempty(R[2][J[2]]) || isempty(R[3][J[3]])
end
function _isempty(R::NTuple{Dim, Vector{UnitRange{Int}}}, J::CartesianIndex{Dim}) where Dim
    return isempty(R[1][J[1]]) || _isempty(R[2:end], CartesianIndex(J.I[2:end]))
end

function refinement_R(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis_R.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges_R.(A, P′)
    return _refinement(M, P′, A, R)
end
function refinement_I(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis_I.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges_I.(A, P′)
    return _refinement(M, P′, A, R)
end
function refinement(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges.(A, P′)
    return _refinement(M, P′, A, R)
end

function _refinement(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}, A::NTuple{Dim, SparseMatrixCSC{U, Int32}}, R::NTuple{Dim, Vector{UnitRange{Int64}}}) where {Dim, Deg, C, T, T′, U}
    a::Array{C, Dim} = controlpoints(M)
    J::CartesianIndex{Dim} = CartesianIndex(findfirst.(!isempty, R))
    D::CartesianIndices{Dim, NTuple{Dim, UnitRange{Int64}}} = CartesianIndices(getindex.(R, J.I))
    value = sum(*(getindex.(A, I.I, J.I)...) * a[I] for I in D)
    j = prod(J.I)
    a′ = __control_points(value, j, dim.(P′))
    for J in view(CartesianIndices(UnitRange.(1, dim.(P′))), (j+1):prod(dim.(P′)))
        if _isempty(R, J)
            @inbounds a′[J] = zero(value)
        else
            D = CartesianIndices(getindex.(R, J.I))
            @inbounds a′[J] = sum(*(getindex.(A, I.I, J.I)...) * a[I] for I in D)
        end
    end
    return BSplineManifold(a′, P′)
end

function refinement_R(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, W, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis_R.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges_R.(A, P′)
    return _refinement(M, P′, A, R)
end
function refinement_I(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, W, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis_I.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges_I.(A, P′)
    return _refinement(M, P′, A, R)
end
function refinement(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}) where {Dim, Deg, C, W, T, T′}
    U = StaticArrays.arithmetic_closure(promote_type(T,T′))
    A::NTuple{Dim, SparseMatrixCSC{U, Int32}} = changebasis.(bsplinespaces(M), P′)
    R::NTuple{Dim, Vector{UnitRange{Int64}}} = _i_ranges.(A, P′)
    return _refinement(M, P′, A, R)
end

function _refinement(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where p′}, A::NTuple{Dim, SparseMatrixCSC{U, Int32}}, R::NTuple{Dim, Vector{UnitRange{Int64}}}) where {Dim, Deg, C, W, T, T′, U}
    a::Array{C, Dim} = controlpoints(M)
    w::Array{W, Dim} = weights(M)
    J::CartesianIndex{Dim} = CartesianIndex(findfirst.(!isempty, R))
    D::CartesianIndices{Dim, NTuple{Dim, UnitRange{Int64}}} = CartesianIndices(getindex.(R, J.I))
    value_w = sum(*(getindex.(A, I.I, J.I)...) * w[I] for I in D)
    value_a = sum(*(getindex.(A, I.I, J.I)...) * w[I] * a[I] for I in D)
    j = prod(J.I)
    w′ = __control_points(value_w, j, dim.(P′))
    a′ = __control_points(value_a, j, dim.(P′))
    for J in view(CartesianIndices(UnitRange.(1, dim.(P′))), (j+1):prod(dim.(P′)))
        if _isempty(R, J)
            @inbounds w′[J] = zero(value_w)
            @inbounds a′[J] = zero(value_a)
        else
            D = CartesianIndices(getindex.(R, J.I))
            @inbounds w′[J] = sum(*(getindex.(A, I.I, J.I)...) * w[I] for I in D)
            @inbounds a′[J] = sum(*(getindex.(A, I.I, J.I)...) * w[I] * a[I] for I in D)
        end
    end
    nans = .!(iszero.(a′) .& iszero.(w′))
    a′ ./= w′
    a′ .*= nans
    return RationalBSplineManifold(a′, w′, P′)
end

function refinement_R(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, T}
    _P′ = _promote_knottype(P′)
    return refinement_R(M, _P′)
end
function refinement_R(M::BSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement_R(M, P′)
end
function refinement_R(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, W, T}
    _P′ = _promote_knottype(P′)
    return refinement_R(M, _P′)
end
function refinement_R(M::RationalBSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement_R(M, P′)
end
function refinement_I(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, T}
    _P′ = _promote_knottype(P′)
    return refinement_I(M, _P′)
end
function refinement_I(M::BSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement_I(M, P′)
end
function refinement_I(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, W, T}
    _P′ = _promote_knottype(P′)
    return refinement_I(M, _P′)
end
function refinement_I(M::RationalBSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement_I(M, P′)
end
function refinement(M::BSplineManifold{Dim, Deg, C, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, T}
    _P′ = _promote_knottype(P′)
    return refinement(M, _P′)
end
function refinement(M::BSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement(M, P′)
end
function refinement(M::RationalBSplineManifold{Dim, Deg, C, W, T}, P′::NTuple{Dim, BSplineSpace{p′,T′} where {p′, T′}}) where {Dim, Deg, C, W, T}
    _P′ = _promote_knottype(P′)
    return refinement(M, _P′)
end
function refinement(M::RationalBSplineManifold{Dim}, P′::Vararg{BSplineSpace, Dim}) where Dim
    return refinement(M, P′)
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knotvector.
"""
@generated function refinement_R(M::AbstractManifold{Dim},
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
    exs = [:($(Symbol(:P,i,"′")) = expandspace_R($(Symbol(:P,i)), $(Symbol(:p,i,:₊)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    return Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exp = p₊),
        :($exk = k₊),
        exs...,
        :(return refinement_R(M, $(exP′)))
    )
end

@generated function refinement_R(M::AbstractManifold{Dim},
                               k₊::NTuple{Dim, AbstractKnotVector}=ntuple(i->EmptyKnotVector(), Val(Dim))) where Dim
    Ps = [Symbol(:P,i) for i in 1:Dim]
    Ps′ = [Symbol(:P,i,"′") for i in 1:Dim]
    ks = [Symbol(:k,i,:₊) for i in 1:Dim]
    exP = Expr(:tuple, Ps...)
    exP′ = Expr(:tuple, Ps′...)
    exk = Expr(:tuple, ks...)
    exs = [:($(Symbol(:P,i,"′")) = expandspace_R($(Symbol(:P,i)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    return Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exk = k₊),
        exs...,
        :(return refinement_R(M, $(exP′)))
    )
end

@doc raw"""
Refinement of B-spline manifold with additional degree and knotvector.
"""
@generated function refinement_I(M::AbstractManifold{Dim},
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
    exs = [:($(Symbol(:P,i,"′")) = expandspace_I($(Symbol(:P,i)), $(Symbol(:p,i,:₊)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    return Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exp = p₊),
        :($exk = k₊),
        exs...,
        :(return refinement_I(M, $(exP′)))
    )
end

@generated function refinement_I(M::AbstractManifold{Dim},
                               k₊::NTuple{Dim, AbstractKnotVector}=ntuple(i->EmptyKnotVector(), Val(Dim))) where Dim
    Ps = [Symbol(:P,i) for i in 1:Dim]
    Ps′ = [Symbol(:P,i,"′") for i in 1:Dim]
    ks = [Symbol(:k,i,:₊) for i in 1:Dim]
    exP = Expr(:tuple, Ps...)
    exP′ = Expr(:tuple, Ps′...)
    exk = Expr(:tuple, ks...)
    exs = [:($(Symbol(:P,i,"′")) = expandspace_I($(Symbol(:P,i)), $(Symbol(:k,i,:₊)))) for i in 1:Dim]
    return Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exk = k₊),
        exs...,
        :(return refinement_I(M, $(exP′)))
    )
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
    return Expr(
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
    return Expr(
        :block,
        :($exP = bsplinespaces(M)),
        :($exk = k₊),
        exs...,
        :(return refinement(M, $(exP′)))
    )
end

# resolve ambiguities
refinement_R(M::AbstractManifold{0}, ::Tuple{}) = M
refinement_R(M::BSplineManifold{0, Deg, C, T, S} where {Deg, C, T, S<:Tuple{}}, ::Tuple{}) = M
refinement_R(M::RationalBSplineManifold{0, Deg, C, W, T, S} where {Deg, C, W, T, S<:Tuple{}}, ::Tuple{}) = M
refinement_I(M::AbstractManifold{0}, ::Tuple{}) = M
refinement_I(M::BSplineManifold{0, Deg, C, T, S} where {Deg, C, T, S<:Tuple{}}, ::Tuple{}) = M
refinement_I(M::RationalBSplineManifold{0, Deg, C, W, T, S} where {Deg, C, W, T, S<:Tuple{}}, ::Tuple{}) = M
refinement(M::AbstractManifold{0}, ::Tuple{}) = M
refinement(M::BSplineManifold{0, Deg, C, T, S} where {Deg, C, T, S<:Tuple{}}, ::Tuple{}) = M
refinement(M::RationalBSplineManifold{0, Deg, C, W, T, S} where {Deg, C, W, T, S<:Tuple{}}, ::Tuple{}) = M

function expand_domain(M::RationalBSplineManifold{2}, Δt::Real)
    if BasicBSpline.isclamped(M)
        M = clamp(M)
    end
    P1,P2 = bsplinespaces(M)
    P1′,P2′ = expand_domain(P1,Δt), expand_domain(P2,Δt)
    P1′′,P2′′ = P1+P1′,P2+P2′
    p1,p2 = degree(P1), degree(P2)
    n1,n2 = dim(P1), dim(P2)
    A1 = inv(Matrix(changebasis_I(P1′, P1′′)[:,begin+p1+1:end-p1-1]))
    A2 = inv(Matrix(changebasis_I(P2′, P2′′)[:,begin+p2+1:end-p2-1]))
    a = controlpoints(M)
    w = weights(M)
    w′ = [sum(A1[i1,j1]*A2[i2,j2]*w[i1,i2] for i1 in 1:n1, i2 in 1:n2) for j1 in 1:n1, j2 in 1:n2]
    a′ = [sum(A1[i1,j1]*A2[i2,j2]*w[i1,i2]*a[i1,i2] for i1 in 1:n1, i2 in 1:n2) for j1 in 1:n1, j2 in 1:n2] ./ w′
    return RationalBSplineManifold(a′,w′,(P1′,P2′))
end

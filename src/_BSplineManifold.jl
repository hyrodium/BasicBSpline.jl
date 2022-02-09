# B-spline manifold

abstract type AbstractBSplineManifold{Dim, Deg} end

dim(::AbstractBSplineManifold{Dim}) where Dim = Dim

@inline function (M::AbstractBSplineManifold{1})(t1)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    unsafe_mapping(M,t1)
end

@inline function (M::AbstractBSplineManifold{2})(t1,t2)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    unsafe_mapping(M,t1,t2)
end

@inline function (M::AbstractBSplineManifold{3})(t1,t2,t3)
    Ps = bsplinespaces(M)
    t1 in domain(Ps[1]) || throw(DomainError(t1, "The input $(t1) is out of range."))
    t2 in domain(Ps[2]) || throw(DomainError(t2, "The input $(t2) is out of range."))
    t3 in domain(Ps[3]) || throw(DomainError(t3, "The input $(t3) is out of range."))
    unsafe_mapping(M,t1,t2,t3)
end

# TODO add mappings higher dimensionnal B-spline manifold with @generated macro

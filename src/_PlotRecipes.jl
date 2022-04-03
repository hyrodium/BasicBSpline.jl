# B-spline space
@recipe function f(P::AbstractBSplineSpace{p}) where p
    # TODO fix number of sampling points
    N = 100(p+1)
    k = knotvector(P)
    ts = Float64[]
    for i in 1:length(k)-1
        append!(ts, range(k[i], k[i+1], length=N))
    end
    unique!(ts)
    n = dim(P)
    tt = repeat([ts;NaN],n)
    bs = [bsplinebasis(P,i,t) for t in ts, i in 1:n]
    bb = vec(vcat(bs,zeros(n)'))
    tt, bb
end

# B-spline Curve
@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{2,<:Real}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2)
end

@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{3,<:Real}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

# Rational B-spline Curve
@recipe function f(M::RationalBSplineManifold{1,Deg,<:StaticVector{2,<:Real}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2)
end

@recipe function f(M::RationalBSplineManifold{1,Deg,<:StaticVector{3,<:Real}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

#=
TODO
* BSplineSurface
* BSplineSolid
* RationalBSplineManifold
=#

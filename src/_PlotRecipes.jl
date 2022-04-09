# B-spline space
@recipe function f(P::AbstractBSplineSpace{p}) where p
    # TODO fix number of sampling points
    N = 100
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

# B-spline curve in 2d
@recipe function f(M::Union{BSplineManifold{1,Deg,<:StaticVector{2,<:Real}}, RationalBSplineManifold{1,Deg,<:StaticVector{2,<:Real}}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :gray
        markershape := :circle
        markercolor := :gray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2)
end

# B-spline curve in 3d
@recipe function f(M::Union{BSplineManifold{1,Deg,<:StaticVector{3,<:Real}}, RationalBSplineManifold{1,Deg,<:StaticVector{3,<:Real}}}) where Deg
    # TODO fix number of sampling points
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=300)
    @series begin
        primary := false
        linecolor := :gray
        markershape := :circle
        markercolor := :gray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

# B-spline surface
@recipe function f(M::Union{BSplineManifold{2,Deg,<:StaticVector{3,<:Real}}, RationalBSplineManifold{2,Deg,<:StaticVector{3,<:Real}}}) where Deg
    # TODO fix number of sampling points
    t1_min, t1_max = extrema(domain(bsplinespaces(M)[1]))
    t2_min, t2_max = extrema(domain(bsplinespaces(M)[2]))
    t1s = range(t1_min, t1_max, length=300)
    t2s = range(t2_min, t2_max, length=300)
    @series begin
        primary := false
        linecolor := :gray
        markershape := :circle
        markercolor := :gray
        seriestype := :scatter
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    @series begin
        primary := false
        linecolor := :gray
        markershape := :circle
        markercolor := :gray
        seriestype := :path
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    @series begin
        primary := false
        linecolor := :gray
        markershape := :circle
        markercolor := :gray
        seriestype := :path
        a = controlpoints(M)
        getindex.(a',1), getindex.(a',2), getindex.(a',3)
    end
    ps = M.(t1s,t2s')
    xs = getindex.(ps,1)
    ys = getindex.(ps,2)
    zs = getindex.(ps,3)
    seriestype := :surface
    xs,ys,zs
end

#=
TODO
* BSplineSolid
=#

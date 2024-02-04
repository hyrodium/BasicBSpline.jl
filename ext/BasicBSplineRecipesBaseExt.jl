module BasicBSplineRecipesBaseExt

using BasicBSpline
using RecipesBase
import BasicBSpline.AbstractFunctionSpace
using StaticArrays

# B-spline space
@recipe function f(P::AbstractFunctionSpace)
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

const _Manifold{Dim1, Dim2} = Union{BSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}, RationalBSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}} where Deg

@recipe function f(M::_Manifold{1, 2})
    # TODO fix number of sampling points
    N = 100
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=N)
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

@recipe function f(M::_Manifold{1, 3})
    # TODO fix number of sampling points
    N = 100
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=N)
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

@recipe function f(M::_Manifold{2, 3})
    # TODO fix number of sampling points
    N = 100
    t1_min, t1_max = extrema(domain(bsplinespaces(M)[1]))
    t2_min, t2_max = extrema(domain(bsplinespaces(M)[2]))
    t1s = range(t1_min, t1_max, length=N)
    t2s = range(t2_min, t2_max, length=N)
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

end # module

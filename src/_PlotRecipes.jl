@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{2,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t_min, t_max, length=300))
    getindex.(p,1), getindex.(p,2)
end

@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{3,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t_min, t_max, length=300))
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

@recipe function f(M::RationalBSplineManifold{1,Deg,<:StaticVector{2,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t_min, t_max, length=300))
    getindex.(p,1), getindex.(p,2)
end

@recipe function f(M::RationalBSplineManifold{1,Deg,<:StaticVector{3,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t_min, t_max, length=300))
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

#=
TODO
* BSplineSurface
* BSplineSolid
* RationalBSplineManifold
=#

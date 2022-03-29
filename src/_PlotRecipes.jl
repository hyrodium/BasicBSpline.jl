@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{2,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1),getindex.(a,2)
    end
    t1,t2 = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t1,t2,length=300))
    getindex.(p,1),getindex.(p,2)
end

@recipe function f(M::BSplineManifold{1,Deg,<:StaticVector{3,<:Real}}) where Deg
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1),getindex.(a,2),getindex.(a,3)
    end
    t1,t2 = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t1,t2,length=300))
    getindex.(p,1),getindex.(p,2),getindex.(p,3)
end

#=
TODO
* BSplineSurface
* BSplineSolid
* RationalBSplineManifold
=#

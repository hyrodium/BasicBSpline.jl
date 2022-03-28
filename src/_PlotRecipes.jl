@recipe function f(M::BSplineManifold{1,ps,<:StaticVector{2,<:Real}}) where ps
    @series begin
        primary := false
        linecolor := :lightgray
        markershape := :circle
        markercolor := :lightgray
        a = controlpoints(M)
        getindex.(a,1),getindex.(a,2)
    end
    t1,t2 = extrema(domain(bsplinespaces(M)[1]))
    p = M.(range(t1,t2,300))
    getindex.(p,1),getindex.(p,2)
end

#=
TODO
* BSplineSurface
* BSplineSolid
* RationalBSplineManifold
=#

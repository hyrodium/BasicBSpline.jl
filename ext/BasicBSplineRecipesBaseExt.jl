module BasicBSplineRecipesBaseExt

using BasicBSpline
using RecipesBase
import BasicBSpline.AbstractFunctionSpace
using StaticArrays

# B-spline space
@recipe function f(P::AbstractFunctionSpace; division_number=10)
    k = knotvector(P)
    ts = Float64[]
    for i in 1:length(k)-1
        append!(ts, range(k[i], k[i+1], length=division_number+1))
    end
    unique!(ts)
    n = dim(P)
    tt = repeat([ts;NaN],n)
    bs = [bsplinebasis(P,i,t) for t in ts, i in 1:n]
    bb = vec(vcat(bs,zeros(n)'))
    tt, bb
end

const _Manifold{Dim1, Dim2} = Union{BSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}, RationalBSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}} where Deg

@kwdef struct PlotAttributesContolPoints
    # https://docs.juliaplots.org/latest/generated/attributes_series/
    line_z=nothing
    linealpha=nothing
    linecolor=:gray
    linestyle=:solid
    linewidth=:auto
    marker_z=nothing
    markeralpha=nothing
    markercolor=:gray
    markershape=:circle
    markersize=4
    markerstrokealpha=nothing
    markerstrokecolor=:match
    markerstrokestyle=:solid
    markerstrokewidth=1
end

@recipe function f(M::_Manifold{1, 2}; controlpoints=(;), division_number=100)
    attributes = PlotAttributesContolPoints(;controlpoints...)
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=division_number+1)
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        marker_z := attributes.marker_z
        markeralpha := attributes.markeralpha
        markercolor := attributes.markercolor
        markershape := attributes.markershape
        markersize := attributes.markersize
        markerstrokealpha := attributes.markerstrokealpha
        markerstrokecolor := attributes.markerstrokecolor
        markerstrokestyle := attributes.markerstrokestyle
        markerstrokewidth := attributes.markerstrokewidth
        a = BasicBSpline.controlpoints(M)
        getindex.(a,1), getindex.(a,2)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2)
end

@recipe function f(M::_Manifold{1, 3}; controlpoints=(;), division_number=100)
    attributes = PlotAttributesContolPoints(;controlpoints...)
    t_min, t_max = extrema(domain(bsplinespaces(M)[1]))
    ts = range(t_min, t_max, length=division_number+1)
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        marker_z := attributes.marker_z
        markeralpha := attributes.markeralpha
        markercolor := attributes.markercolor
        markershape := attributes.markershape
        markersize := attributes.markersize
        markerstrokealpha := attributes.markerstrokealpha
        markerstrokecolor := attributes.markerstrokecolor
        markerstrokestyle := attributes.markerstrokestyle
        markerstrokewidth := attributes.markerstrokewidth
        a = BasicBSpline.controlpoints(M)
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    p = M.(ts)
    getindex.(p,1), getindex.(p,2), getindex.(p,3)
end

@recipe function f(M::_Manifold{2, 2}; controlpoints=(;), division_number=100)
    attributes = PlotAttributesContolPoints(;controlpoints...)
    a = BasicBSpline.controlpoints(M)
    @series begin
        primary := false
        marker_z := attributes.marker_z
        markeralpha := attributes.markeralpha
        markercolor := attributes.markercolor
        markershape := attributes.markershape
        markersize := attributes.markersize
        markerstrokealpha := attributes.markerstrokealpha
        markerstrokecolor := attributes.markerstrokecolor
        markerstrokestyle := attributes.markerstrokestyle
        markerstrokewidth := attributes.markerstrokewidth
        seriestype := :scatter
        getindex.(a,1), getindex.(a,2)
    end
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        seriestype := :path
        getindex.(a,1), getindex.(a,2)
    end
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        seriestype := :path
        getindex.(a',1), getindex.(a',2)
    end
    I1, I2 = domain.(bsplinespaces(M))
    a1, b1 = float.(extrema(I1))
    a2, b2 = float.(extrema(I2))
    ts_boundary = [
        [(t1, a2) for t1 in view(range(I1, length=division_number+1), 1:division_number)];
        [(b1, t2) for t2 in view(range(I2, length=division_number+1), 1:division_number)];
        [(t1, b2) for t1 in view(range(I1, length=division_number+1), division_number+1:-1:2)];
        [(a1, t2) for t2 in view(range(I2, length=division_number+1), division_number+1:-1:1)]
    ]
    ps_boundary = [M(t...) for t in ts_boundary]
    fill := true
    fillalpha --> 0.5
    getindex.(ps_boundary,1), getindex.(ps_boundary,2)
end

@recipe function f(M::_Manifold{2, 3}; controlpoints=(;), division_number=100)
    attributes = PlotAttributesContolPoints(;controlpoints...)
    t1_min, t1_max = extrema(domain(bsplinespaces(M)[1]))
    t2_min, t2_max = extrema(domain(bsplinespaces(M)[2]))
    t1s = range(t1_min, t1_max, length=division_number+1)
    t2s = range(t2_min, t2_max, length=division_number+1)
    a = BasicBSpline.controlpoints(M)
    @series begin
        primary := false
        marker_z := attributes.marker_z
        markeralpha := attributes.markeralpha
        markercolor := attributes.markercolor
        markershape := attributes.markershape
        markersize := attributes.markersize
        markerstrokealpha := attributes.markerstrokealpha
        markerstrokecolor := attributes.markerstrokecolor
        markerstrokestyle := attributes.markerstrokestyle
        markerstrokewidth := attributes.markerstrokewidth
        seriestype := :scatter
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        seriestype := :path
        getindex.(a,1), getindex.(a,2), getindex.(a,3)
    end
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        seriestype := :path
        getindex.(a',1), getindex.(a',2), getindex.(a',3)
    end
    ps = M.(t1s,t2s')
    xs = getindex.(ps,1)
    ys = getindex.(ps,2)
    zs = getindex.(ps,3)
    seriestype := :surface
    xs, ys, zs
end

#=
TODO
* BSplineSolid
=#

end # module

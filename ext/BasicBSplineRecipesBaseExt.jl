module BasicBSplineRecipesBaseExt

using BasicBSpline
using RecipesBase
import BasicBSpline.AbstractFunctionSpace
using StaticArrays

# B-spline space
@recipe function f(k::AbstractKnotVector; shift_y=0.02, plane=:xy)
    u = BasicBSpline._vec(k)
    l = length(u)
    v = zeros(Float64, l)
    for i in 2:l
        if u[i-1] == u[i]
            v[i] = v[i-1] - shift_y
        end
    end
    seriestype := :scatter
    delete!(plotattributes, :shift_y)
    delete!(plotattributes, :plane)
    if plane === :xy
        u, v
    elseif plane === :xy
        v, u
    elseif plane === :xz
        u, zero(u), v
    elseif plane === :yz
        zero(u), u, v
    end
end

# B-spline space
@recipe function f(P::AbstractFunctionSpace; division_number=20, plane=:xy)
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
    delete!(plotattributes, :division_number)
    delete!(plotattributes, :plane)
    if plane === :xy
        tt, bb
    elseif plane === :xy
        bb, tt
    elseif plane === :xz
        tt, zero(tt), bb
    elseif plane === :yz
        zero(tt), tt, bb
    end
end

const _Manifold{Dim1, Dim2} = Union{BSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}, RationalBSplineManifold{Dim1,Deg,<:StaticVector{Dim2,<:Real}}} where Deg

Base.@kwdef struct PlotAttributesContolPoints
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
    delete!(plotattributes, :controlpoints)
    delete!(plotattributes, :division_number)
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
    delete!(plotattributes, :controlpoints)
    delete!(plotattributes, :division_number)
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
    delete!(plotattributes, :controlpoints)
    delete!(plotattributes, :division_number)
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
    delete!(plotattributes, :controlpoints)
    delete!(plotattributes, :division_number)
    xs, ys, zs
end

@recipe function f(M::_Manifold{3, 3}; controlpoints=(;), division_number=100)
    attributes = PlotAttributesContolPoints(;controlpoints...)
    t1_min, t1_max = extrema(domain(bsplinespaces(M)[1]))
    t2_min, t2_max = extrema(domain(bsplinespaces(M)[2]))
    t3_min, t3_max = extrema(domain(bsplinespaces(M)[3]))
    t1s = range(t1_min, t1_max, length=division_number+1)
    t2s = range(t2_min, t2_max, length=division_number+1)
    t3s = range(t3_min, t3_max, length=division_number+1)
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
        vec(getindex.(a,1)), vec(getindex.(a,2)), vec(getindex.(a,3))
    end
    @series begin
        primary := false
        line_z := attributes.line_z
        linealpha := attributes.linealpha
        linecolor := attributes.linecolor
        linestyle := attributes.linestyle
        linewidth := attributes.linewidth
        seriestype := :path
        n1, n2, n3 = dim.(bsplinespaces(M))
        X = Float64[]
        Y = Float64[]
        Z = Float64[]
        a1 = a
        xs = [p[1] for p in a1]
        ys = [p[2] for p in a1]
        zs = [p[3] for p in a1]
        append!(X, vec(cat(xs, fill(NaN,1,n2,n3), dims=1)))
        append!(Y, vec(cat(ys, fill(NaN,1,n2,n3), dims=1)))
        append!(Z, vec(cat(zs, fill(NaN,1,n2,n3), dims=1)))
        a2 = permutedims(a,(2,3,1))
        xs = [p[1] for p in a2]
        ys = [p[2] for p in a2]
        zs = [p[3] for p in a2]
        append!(X, vec(cat(xs, fill(NaN,1,n3,n1), dims=1)))
        append!(Y, vec(cat(ys, fill(NaN,1,n3,n1), dims=1)))
        append!(Z, vec(cat(zs, fill(NaN,1,n3,n1), dims=1)))
        a3 = permutedims(a,(3,1,2))
        xs = [p[1] for p in a3]
        ys = [p[2] for p in a3]
        zs = [p[3] for p in a3]
        append!(X, vec(cat(xs, fill(NaN,1,n1,n2), dims=1)))
        append!(Y, vec(cat(ys, fill(NaN,1,n1,n2), dims=1)))
        append!(Z, vec(cat(zs, fill(NaN,1,n1,n2), dims=1)))

        X,Y,Z
    end

    nanvec = fill(SVector(NaN,NaN,NaN), 1, division_number+1)

    qs1 = M.(t1s, t2s', t3_max)
    qs2 = M.(t1s', t2s, t3_min)
    qs3 = M.(t1_max, t2s, t3s')
    qs4 = M.(t1_min, t2s', t3s)
    qs5 = M.(t1s', t2_max, t3s)
    qs6 = M.(t1s, t2_min, t3s')

    ps = [
        qs1;
        qs1[end:end,:];  # Adhoc additional values for fixing boundary
        nanvec;
        qs2;
        qs2[end:end,:];
        nanvec;
        qs3;
        qs3[end:end,:];
        nanvec;
        qs4;
        qs4[end:end,:];
        nanvec;
        qs5;
        qs5[end:end,:];
        nanvec;
        qs6
    ]

    xs = getindex.(ps,1)
    ys = getindex.(ps,2)
    zs = getindex.(ps,3)

    seriestype := :surface
    delete!(plotattributes, :controlpoints)
    delete!(plotattributes, :division_number)
    xs, ys, zs
end

end # module

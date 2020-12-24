# Geometry

"""
tangentvectors on B-spline manifold
NOTICE: these methods are experimental
"""
tangentvectors

function tangentvectors(M::BSplineCurve{p1}, t1::Real) where {p1}
    P1, = bsplinespaces(M)
    a = controlpoints(M)
    n1 = dim(P1)
    k1 = P1.knots
    j1 = _knotindex(P1, t1)
    b1′ = _bsplinebasis′(P1,t1,j1)
    d̂ = size(a)[end]

    S1 = Array{Float64}(undef, d̂)
    f1 = j1-p1
    l1 = j1
    for xx in 1:d̂
        S1[xx] = b1′[1]*a[f1,xx]
        for i1 in f1+1:l1
            S1[xx] += b1′[i1-f1+1]*a[i1,xx]
        end
    end

    return S1
end


function tangentvectors(M::BSplineSurface{p1,p2}, t1::Real,t2::Real) where {p1} where {p2}
    P1, P2 = bsplinespaces(M)
    a = controlpoints(M)
    n1, n2 = dim(P1), dim(P2)
    k1, k2 = P1.knots, P2.knots
    j1, j2 = _knotindex(P1, t1), _knotindex(P2, t2)
    b1 = _bsplinebasis(P1,t1,j1)
    b2 = _bsplinebasis(P2,t2,j2)
    b1′ = _bsplinebasis′(P1,t1,j1)
    b2′ = _bsplinebasis′(P2,t2,j2)
    d̂ = size(a)[end]

    S1 = Array{Float64}(undef, d̂)
    S2 = Array{Float64}(undef, d̂)
    f1, f2 = j1-p1, j2-p2
    l1, l2 = j1, j2

    for xx in 1:d̂
        S2[xx] = b2[1]*a[f1,f2,xx]
        for i2 in f2+1:l2
            S2[xx] += b2[i2-f2+1]*a[f1,i2,xx]
        end
        S1[xx] = b1′[1]*S2[xx]
        for i1 in f1+1:l1
            S2[xx] = b2[1]*a[i1,f2,xx]
            for i2 in f2+1:l2
                S2[xx] += b2[i2-f2+1]*a[i1,i2,xx]
            end
            S1[xx] += b1′[i1-f1+1]*S2[xx]
        end
    end
    X1 = copy(S1)

    for xx in 1:d̂
        S2[xx] = b2′[1]*a[f1,f2,xx]
        for i2 in f2+1:l2
            S2[xx] += b2′[i2-f2+1]*a[f1,i2,xx]
        end
        S1[xx] = b1[1]*S2[xx]
        for i1 in f1+1:l1
            S2[xx] = b2′[1]*a[i1,f2,xx]
            for i2 in f2+1:l2
                S2[xx] += b2′[i2-f2+1]*a[i1,i2,xx]
            end
            S1[xx] += b1[i1-f1+1]*S2[xx]
        end
    end
    X2 = copy(S1)

    return (X1,X2)
end

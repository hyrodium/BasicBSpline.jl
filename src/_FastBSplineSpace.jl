# B-spline space
struct FastBSplineSpace{T}
    vector::Array{Float64,1}
    function FastBSplineSpace(p::Int, knots::Knots)
        if p < 0
            error("degree of polynominal must be non-negative")
        elseif p > MAX_DEGREE
            error("FastBSpline supports only degree 0 , ... , 3")
        end
        new{p}(knots.vector)
    end
end

function Knots(P::FastBSplineSpace)
    return Knots(P.vector)
end


for p in 0:MAX_DEGREE
    @eval begin
        @doc raw"""
        Return dimention of a B-spline space.
        ```math
        \dim(\mathcal{P}[p,k])
        =\sharp k - p -1
        ```
        """
        function dim(P::FastBSplineSpace{$p})
            k=P.vector
            return length(k)-$p-1
        end
    end

    @eval begin
        @doc raw"""
        Return dimention of a B-spline space.
        ```math
        \dim(\mathcal{P}[p,k])
        =\sharp k - p -1
        ```
        """
        function BSplineSpace(P::FastBSplineSpace{$p})
            BSplineSpace($p, Knots(P.vector))
        end
    end
end


for p in 0:MAX_DEGREE, p′ in 0:MAX_DEGREE
    @eval begin
        """
        Check inclusive relationship between B-spline spaces.
        """
        function Base.:⊆(P::FastBSplineSpace{$p}, P′::FastBSplineSpace{$p′})
            k = Knots(P)
            k′ = Knots(P′)
            p₊ = $p′-$p

            return (k+p₊*unique(k) ⊆ k′) && p₊ ≥ 0
        end
    end
end


for p in 0:MAX_DEGREE
    @eval begin
        function iszeros(P::FastBSplineSpace{$p})
            k = P.vector
            n = dim(P)
            return [k[i] == k[i+$p+1] for i ∈ 1:n]
        end
    end
end

function isproper(P::FastBSplineSpace)
    return !|(iszeros(P)...)
end

function properdim(P::FastBSplineSpace)
    return dim(P) - sum(iszeros(P))
end

@testset "compatibility" begin
    p = 2
    k = Knots(rand(10)) + p * Knots(0,1)
    P = FastBSplineSpace(p, k)
    n = dim(P)
    i = n÷2
    t = 0.5

    bsplinebasis₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinebasis₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinebasis(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinebasis′₊₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinebasis′₋₀(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinebasis′(i::Integer,P::AbstractBSplineSpace,t::Real) = bsplinebasis(P,i,t)
    bsplinesupport(i::Integer,P::AbstractBSplineSpace) = bsplinesupport(P,i)
end

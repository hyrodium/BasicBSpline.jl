using BenchmarkTools
using BasicBSpline
using GeometryBasics
using Random

const SUITE = BenchmarkGroup()

Random.seed!(42)

p = 3
k = Knots(rand(100)) + (p+1) * Knots(0,1)
P = FastBSplineSpace(p,k)
i = 1
t = 0.02

SUITE["fast-bsplinebasis"] = @benchmarkable bsplinebasis(i,P,t)

p1 = 2
k1 = Knots(rand(3)) + (p1 + 1) * Knots([0, 1])
P1 = FastBSplineSpace(p1, k1)
n1 = dim(P1)
a = [Point(i1, rand()) for i1 in 1:n1]
M = BSplineCurve([P1], a)

p1′ = p1 + 1
k1′ = k1 + unique(k1) + Knots(rand(2))
P1′ = FastBSplineSpace(p1′, k1′)
M′ = refinement(M, [P1′])

SUITE["refinement"] = @benchmarkable refinement(M, [P1′])
SUITE["fitting_I"] = @benchmarkable fittingcontrolpoints(M, [P1′])
SUITE["fitting_R"] = @benchmarkable fittingcontrolpoints(M, [P1′], domain=:R)

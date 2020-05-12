using BasicBSpline
using BenchmarkTools

k=Knots(rand(8))

P=BSplineSpace(3,k)
fP=FastBSplineSpace(3,k)

dim(P)
dim(fP)

bsplinebasis(1,P,0.2)
bsplinebasis(1,fP,0.2)

@benchmark bsplinebasis(1,P,0.4)
@benchmark bsplinebasis(1,fP,0.4)

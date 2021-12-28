using BasicBSpline
using Plots
gr()

##
p = 2
k = KnotVector(1:8)
P = BSplineSpace(p, k)
plot([t -> bsplinebasis₊₀(P, i, t) for i in 1:dim(P)], 1, 8, ylims = (0, 1))
savefig("docs/src/img/bsplinebasisplot.png")

##
p = 2
k = KnotVector(1:8)
P = BSplineSpace(p, k)
plot([t -> bsplinebasis′₊₀(P, i, t) for i in 1:dim(P)], 1, 8)
savefig("docs/src/img/bsplinebasisderivativeplot.png")

##
p = 2
k = KnotVector(1:8)
P = BSplineSpace(p, k)
plot(t -> sum(bsplinebasis₊₀(P, i, t) for i in 1:dim(P)), 1, 8, ylims = (0, 1))
savefig("docs/src/img/sumofbsplineplot.png")

##
p = 2
k = KnotVector(1:8) + p * KnotVector([1, 8])
P = BSplineSpace(p, k)
plot(t -> sum(bsplinebasis₊₀(P, i, t) for i in 1:dim(P)), 1, 8, ylims = (0, 1))
savefig("docs/src/img/sumofbsplineplot2.png")

##
p = 2
k = KnotVector(1:8) + p * KnotVector([1, 8])
P = BSplineSpace(p, k)
plot(t -> sum(bsplinebasis(P, i, t) for i in 1:dim(P)), 1, 8, ylims = (0, 1))
savefig("docs/src/img/sumofbsplineplot3.png")

##
P1 = BSplineSpace(1, KnotVector([1, 3, 5, 8]))
P2 = BSplineSpace(1, KnotVector([1, 3, 5, 6, 8, 9]))
P3 = BSplineSpace(2, KnotVector([1, 1, 3, 3, 5, 5, 8, 8]))
plot([t -> bsplinebasis₊₀(P1, i, t) for i in 1:dim(P1)], 1, 9, ylims = (0, 1))
savefig("docs/src/img/bsplineplotP1.png")
plot([t -> bsplinebasis₊₀(P2, i, t) for i in 1:dim(P2)], 1, 9, ylims = (0, 1))
savefig("docs/src/img/bsplineplotP2.png")
plot([t -> bsplinebasis₊₀(P3, i, t) for i in 1:dim(P3)], 1, 9, ylims = (0, 1))
savefig("docs/src/img/bsplineplotP3.png")

plot(
    plot([t -> bsplinebasis₊₀(P1, i, t) for i in 1:dim(P1)], 1, 9, ylims = (0, 1), legend = false),
    plot([t -> bsplinebasis₊₀(P2, i, t) for i in 1:dim(P2)], 1, 9, ylims = (0, 1), legend = false),
    plot([t -> bsplinebasis₊₀(P3, i, t) for i in 1:dim(P3)], 1, 9, ylims = (0, 1), legend = false),
    layout = (3, 1),
    link = :x,
)
savefig("docs/src/img/subbsplineplot.png")

##
A12 = changebasis(P1, P2)
A13 = changebasis(P1, P3)

plot(
    plot([t -> bsplinebasis₊₀(P1, i, t) for i in 1:dim(P1)], 1, 9, ylims = (0, 1), legend = false),
    plot([t -> sum(A12[i, j] * bsplinebasis₊₀(P2, j, t) for j in 1:dim(P2)) for i in 1:dim(P1)], 1, 9, ylims = (0, 1), legend = false),
    plot([t -> sum(A13[i, j] * bsplinebasis₊₀(P3, j, t) for j in 1:dim(P3)) for i in 1:dim(P1)], 1, 9, ylims = (0, 1), legend = false),
    layout = (3, 1),
    link = :x,
)
savefig("docs/src/img/subbsplineplot2.png")

##
p = 3
k = KnotVector([0.00, 1.50, 2.50, 5.50, 8.00, 9.00, 9.50, 10.0])
P0 = FastBSplineSpace(0, k)
P1 = FastBSplineSpace(1, k)
P2 = FastBSplineSpace(2, k)
P3 = FastBSplineSpace(3, k)
plot(
    plot([t -> bsplinebasis(P0, i, t) for i in 1:dim(P0)], 0, 10, ylims = (0, 1), legend = false),
    plot([t -> bsplinebasis(P1, i, t) for i in 1:dim(P1)], 0, 10, ylims = (0, 1), legend = false),
    plot([t -> bsplinebasis(P2, i, t) for i in 1:dim(P2)], 0, 10, ylims = (0, 1), legend = false),
    plot([t -> bsplinebasis(P3, i, t) for i in 1:dim(P3)], 0, 10, ylims = (0, 1), legend = false),
    layout = (2, 2),
)
savefig("docs/src/img/cover.png")

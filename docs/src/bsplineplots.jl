using BasicBSpline
using Plots
gr()

##
p = 2
k = Knots(1:8)
P = BSplineSpace(p,k)
plot([t->bsplinebasis₊₀(i,P,t) for i in 1:dim(P)], 1, 8, ylims=(0,1.05))
savefig("docs/src/img/bsplinebasisplot.png")

##
p = 2
k = Knots(1:8)
P = BSplineSpace(p,k)
plot([t->bsplinebasis′₊₀(i,P,t) for i in 1:dim(P)], 1, 8)
savefig("docs/src/img/bsplinebasisderivativeplot.png")

##
p = 2
k = Knots(1:8)
P = BSplineSpace(p,k)
plot(t->sum(bsplinebasis₊₀(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))
savefig("docs/src/img/sumofbsplineplot.png")

##
p = 2
k = Knots(1:8) + p * Knots([1,8])
P = BSplineSpace(p,k)
plot(t->sum(bsplinebasis₊₀(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))
savefig("docs/src/img/sumofbsplineplot2.png")

##
p = 2
k = Knots(1:8) + p * Knots([1,8])
P = BSplineSpace(p,k)
plot(t->sum(bsplinebasis(i,P,t) for i in 1:dim(P)), 1, 8, ylims=(0,1.05))
savefig("docs/src/img/sumofbsplineplot3.png")

##
P1 = BSplineSpace(1,Knots([1,3,5,8]))
P2 = BSplineSpace(1,Knots([1,3,5,6,8,9]))
P3 = BSplineSpace(2,Knots([1,1,3,3,5,5,8,8]))
plot([t->bsplinebasis₊₀(i,P1,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05))
savefig("docs/src/img/bsplineplotP1.png")
plot([t->bsplinebasis₊₀(i,P2,t) for i in 1:dim(P2)], 1, 9, ylims=(0,1.05))
savefig("docs/src/img/bsplineplotP2.png")
plot([t->bsplinebasis₊₀(i,P3,t) for i in 1:dim(P3)], 1, 9, ylims=(0,1.05))
savefig("docs/src/img/bsplineplotP3.png")

plot(
    plot([t->bsplinebasis₊₀(i,P1,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),
    plot([t->bsplinebasis₊₀(i,P2,t) for i in 1:dim(P2)], 1, 9, ylims=(0,1.05), legend=false),
    plot([t->bsplinebasis₊₀(i,P3,t) for i in 1:dim(P3)], 1, 9, ylims=(0,1.05), legend=false),
    layout=(3,1),
    link=:x
)
savefig("docs/src/img/subbsplineplot.png")

##
A12 = BasicBSpline.changebasis(P1,P2)
A13 = BasicBSpline.changebasis(P1,P3)

plot(
    plot([t->bsplinebasis₊₀(i,P1,t) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),
    plot([t->sum(A12[i,j]*bsplinebasis₊₀(j,P2,t) for j in 1:dim(P2)) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),
    plot([t->sum(A13[i,j]*bsplinebasis₊₀(j,P3,t) for j in 1:dim(P3)) for i in 1:dim(P1)], 1, 9, ylims=(0,1.05), legend=false),
    layout=(3,1),
    link=:x
)
savefig("docs/src/img/subbsplineplot2.png")

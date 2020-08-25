using BasicBSpline
using ExportNURBS

##

p = 2 # degree of polynomial
k = Knots(1:8) # knot vector
P = BSplineSpace(p,k) # B-spline space
rand_a = [rand(2) for i in 1:dim(P), j in 1:dim(P)]
a = [[2*i-6.5,2*j-6.5] for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
M = BSplineManifold([P,P],a) # Define B-spline manifold
save_png("docs/src/img/2dim.png", M, unitlength=50)

##

k₊=[Knots(3.3,4.2),Knots(3.8,3.2,5.3)]
M′ = refinement(M,k₊=k₊)
save_png("docs/src/img/2dim_refinement.png", M′, unitlength=50)

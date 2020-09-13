using BasicBSpline
using ExportNURBS

## 2-dim B-spline manifold
p = 2 # degree of polynomial
k = Knots(1:8) # knot vector
P = BSplineSpace(p, k) # B-spline space
rand_a = [rand(2) for i in 1:dim(P), j in 1:dim(P)]
a = [[2 * i - 6.5, 2 * j - 6.5] for i in 1:dim(P), j in 1:dim(P)] + rand_a # random generated control points
M = BSplineManifold([P, P], a) # Define B-spline manifold
save_png("docs/src/img/2dim.png", M, unitlength = 50)

## Refinement
k₊ = [Knots(3.3, 4.2), Knots(3.8, 3.2, 5.3)]
M′ = refinement(M, k₊ = k₊)
save_png("docs/src/img/2dim_refinement.png", M′, unitlength = 50)

## Fitting
p1 = 2
p2 = 2
k1 = Knots(-10:10) + p1 * Knots(-10, 10)
k2 = Knots(-10:10) + p2 * Knots(-10, 10)
P1 = FastBSplineSpace(p1, k1)
P2 = FastBSplineSpace(p2, k2)

f(u) = [2u[1] + sin(u[1]) + cos(u[2]) + u[2] / 2, 3u[2] + sin(u[2]) + sin(u[1]) / 2 + u[1]^2 / 6] / 5

a = fittingcontrolpoints(f, [P1, P2])
M = BSplineManifold([P1, P2], a)
save_png("docs/src/img/fitting.png", M, unitlength = 50, up = 10, down = -10, left = -10, right = 10)

## Coarse fitting
p1 = 2
p2 = 2
k1 = Knots(-10:5:10) + p1 * Knots(-10, 10)
k2 = Knots(-10:5:10) + p2 * Knots(-10, 10)
P1 = FastBSplineSpace(p1, k1)
P2 = FastBSplineSpace(p2, k2)

f(u) = [2u[1] + sin(u[1]) + cos(u[2]) + u[2] / 2, 3u[2] + sin(u[2]) + sin(u[1]) / 2 + u[1]^2 / 6] / 5

a = fittingcontrolpoints(f, [P1, P2])
M = BSplineManifold([P1, P2], a)
save_png("docs/src/img/fitting_coarse.png", M, unitlength = 50, up = 10, down = -10, left = -10, right = 10)

## Sine curve
p = 3
k = Knots(range(-2π, 2π, length = 8)) + p * Knots(-2π, 2π)
P = FastBSplineSpace(p, k)

f(u) = [u[1], sin(u[1])]

a = fittingcontrolpoints(f, [P])
M = BSplineManifold([P], a)
save_svg("docs/src/img/sine_curve.svg", M, unitlength = 50, up = 2, down = -2, left = -8, right = 8)

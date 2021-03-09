@testset begin
    gl1 = BasicBSpline.GaussLegendre(1)
    gl2 = BasicBSpline.GaussLegendre(2)
    gl3 = BasicBSpline.GaussLegendre(3)
    gl4 = BasicBSpline.GaussLegendre(4)

    I = 2..8
    coef = rand(10)

    f(n,x) = sum(coef[i]*x^(i-1) for i in 1:n)
    F(n,x) = sum(coef[i]*x^(i)/i for i in 1:n)

    @test F(1,maximum(I)) - F(1,minimum(I)) ≈ BasicBSpline.integrate(x->f(1,x), I, gl1)
    @test F(2,maximum(I)) - F(2,minimum(I)) ≈ BasicBSpline.integrate(x->f(2,x), I, gl1)
    @test F(3,maximum(I)) - F(3,minimum(I)) ≈ BasicBSpline.integrate(x->f(3,x), I, gl2)
    @test F(4,maximum(I)) - F(4,minimum(I)) ≈ BasicBSpline.integrate(x->f(4,x), I, gl2)
    @test F(5,maximum(I)) - F(5,minimum(I)) ≈ BasicBSpline.integrate(x->f(5,x), I, gl3)
    @test F(6,maximum(I)) - F(6,minimum(I)) ≈ BasicBSpline.integrate(x->f(6,x), I, gl3)
    @test F(7,maximum(I)) - F(7,minimum(I)) ≈ BasicBSpline.integrate(x->f(7,x), I, gl4)
    @test F(8,maximum(I)) - F(8,minimum(I)) ≈ BasicBSpline.integrate(x->f(8,x), I, gl4)
end

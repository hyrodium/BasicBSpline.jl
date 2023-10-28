@setup_workload begin
    k1 = KnotVector(1:8)
    k2 = UniformKnotVector(1:8)
    P1 = BSplineSpace{3}(k1)
    P2 = BSplineSpace{3}(k2)
    @compile_workload begin
        changebasis_I(P1, P1)
        changebasis_R(P1, P1)
        changebasis_I(P1, P2)
        changebasis_R(P1, P2)
        changebasis_I(P2, P2)
        changebasis_R(P2, P2)
    end
end

using BasicBSpline
using Plots
using IntervalSets
using Random

function vis_R(P1, P2)
    @assert P1 ⊆ P2
    k1 = knotvector(P1)
    k2 = knotvector(P2)    
    xmin = min(k1[1], k2[1])
    xmax = max(k1[end], k2[end])
    ps = Plots.Plot[]
    A = changebasis_R(P1, P2)
    for i in 1:dim(P1)
        p = plot(t -> bsplinebasis(P1, i, t), extrema(bsplinesupport(P1,i))..., xlims=(xmin, xmax), ylims=(0,1), legend=false, color=:red)
        j_begin, j_end = BasicBSpline._find_j_begin_end_R(P1, P2, i, 1, 1)
        for j in j_begin:j_end
            plot!(p, t -> A[i,j]*bsplinebasis(P2, j, t), extrema(bsplinesupport(P1,i))..., color=:blue)
        end
        push!(ps, p)
    end
    plot(ps..., layout = (dim(P1), 1), size=(700,1000))
end

function vis_I(P1, P2)
    @assert P1 ⊑ P2
    k1 = knotvector(P1)
    k2 = knotvector(P2)    
    xmin = min(k1[1], k2[1])
    xmax = max(k1[end], k2[end])
    ps = Plots.Plot[]
    A = changebasis_I(P1, P2)
    for i in 1:dim(P1)
        p = plot(t -> bsplinebasis(P1, i, t), extrema(bsplinesupport(P1,i))..., xlims=(xmin, xmax), ylims=(-0.5,1), legend=false, color=:red)
        j_begin, j_end = BasicBSpline._find_j_begin_end_I(P1, P2, i, 1, 1)
        for j in j_begin:j_end
            plot!(p, t -> A[i,j]*bsplinebasis(P2, j, t), extrema(bsplinesupport(P2,j))..., color=:blue)
        end
        plot!(p, [extrema(domain(P1))...], seriestype=:vline)
        push!(ps, p)
    end
    plot(ps..., layout = (dim(P1), 1), size=(700,1000))
end

# Rの可視化は簡単だね
Random.seed!(4)
k1 = KnotVector(rand(5)) + 4*KnotVector([0,1])
k2 = k1 + KnotVector(rand(3))
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
vis_R(P1,P2)

# 赤がはみ出したら全部出てくる、が基本
k1 = knotvector" 1 1 11  1 1 1 1"
k2 = knotvector"1 1 1 1  11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 赤がはみ出したら全部出てくる、が基本 (2)
k1 = knotvector" 1 1 11   1 1 1 1"
k2 = knotvector"1 1 1 1 2 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 赤がはみ出ても、多少一致してたら少なくなる
k1 = knotvector"11 11   1 1 1 1"
k2 = knotvector"11 11 2 11   1  1"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 赤がはみ出ても、多少一致してたら少なくなる (2)
k1 = knotvector"11 11   22"
k2 = knotvector"11 11 2 22"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# もう少し詳しく。内側から1つずれた場合 → 全部でてくる (下三角 3)
k1 = knotvector"11 11   22"
k2 = knotvector"111 1  222"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 内側の2つ目がずれた場合 → 多少でてくる (下三角 2)
k1 = knotvector"11 11   22"
k2 = knotvector"1 111  222"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 逆向きに内側の2つ目がずれた場合 → 多少でてくる (下三角 2)
k1 = knotvector"1 111   22"
k2 = knotvector"11 11  222"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
changebasis(P1,P2)
vis_I(P1,P2)

# 予想: 多項式次数が等しい場合、k[2]...k[p]のズレによってブロックの数がずれる
# 多項式次数が等しくない場合、縮退がある場合、混迷を極める…

# でっけぇ空間で包み込み作戦も失敗…
k1 = knotvector"1 111   22"
k2 = knotvector"11 11  222"
k3 = knotvector"11111  222"
P1 = BSplineSpace{3}(k1)
P2 = BSplineSpace{3}(k2)
P3 = BSplineSpace{3}(k3)

P1 ⊑ P2
P2 ⊆ P3
P1 ⊆ P3

A = BasicBSpline._changebasis_I(P1, P2)
B = BasicBSpline._changebasis_R(P2, P3)
C = BasicBSpline._changebasis_R(P1, P3)

A*B - C

(!iszero).(A) * (!iszero).(B)
(!iszero).(C)

# ではどうするか?
# ノット列の両端のズレは諦める
#   そもそもissue解決できてないやん
# 頑張って法則を見つける
#   もういい加減mergeしたい
# ノット列を単純化して、左端と右端の始末をする
#   愚直だが最悪ではない。(浮動小数点の安定性に関して)
#   NaNの問題って解決できるんけ?



BasicBSpline._changebasis_I(P1, P2)

Basic


using BasicBSpline
Pd1 = BSplineSpace{3}(knotvector"4  121 ")
Pd2 = BSplineSpace{4}(knotvector"51 2121")
BasicBSpline._changebasis_I_old(Pd1, Pd2)
Matrix(BasicBSpline._changebasis_I_old(Pd1, Pd2))

Pd1 = BSplineSpace{2}(knotvector"3  12 ")
Pd2 = BSplineSpace{3}(knotvector"4  212")
BasicBSpline._changebasis_I_old(Pd1, Pd2)


using BasicBSpline
Pd1 = BSplineSpace{2}(knotvector" 21 3")
Pd2 = BSplineSpace{3}(knotvector"212 4")
BasicBSpline.changebasis_I(Pd1, Pd2)
isdegenerate_R(Pd1)
isdegenerate_R(Pd2)
isdegenerate_I(Pd1)
isdegenerate_I(Pd2)



s = " 1 111 23"
sum(KnotVector(findall(==('0'+i), s))*i for i in 1:9)

testcases = Tuple{BSplineSpace, BSplineSpace}[]
p = 3
chars = '0' .+ (0:p+1)
for _ in 1:1000000
    s1 = join(rand(chars,5))
    s2 = join(rand(chars,5))
    k1 = sum(KnotVector(findall(==('0'+i), s1))*i for i in 1:9)
    k2 = sum(KnotVector(findall(==('0'+i), s2))*i for i in 1:9)
    P1 = BSplineSpace{p}(k1)
    P2 = BSplineSpace{p}(k2)
    dim(P1) < 1 && continue
    dim(P2) < 1 && continue
    P1 ⊑ P2 && push!(testcases, (P1, P2))
    P2 ⊑ P1 && push!(testcases, (P2, P1))
end

testcases

unique(testcases)

P1, P2 = testcases[1]

dim(P2)

plot(P1)
plot!(P2)

changebasis(P1, P2)
BasicBSpline.changebasis_I(P1, P2)
BasicBSpline._changebasis_I_old(P1, P2)

dim(BSplineSpace{3}(KnotVector([4, 4, 4])))

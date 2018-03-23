include("NMF.jl")

using NMF

# 通常版を使用する場合

# 分解対象の[2,4]の行列を用意
V = [1 2 3 4; 2 3 4 5]

# 因子分回数 K を５で設定する
k = 3

# NMF処理を実行
# @time W,H = nmf_euc(V,k)
# @time W,H = nmf_kl(V,k)
@time W,H = nmf_is(V,k)

# 結果を確認
println("Normal NMF Result : $(norm(V - W*H))")
println("-------------------------------------------")


# テンプレート版を使用する場合

# 分解対象の[2,4]の行列を用意
V = [1 2 3 4; 2 3 4 5]

# テンプレート行列を用意
W = [0. 1. 11.;2. 6. 10.]

# NMF処理を実行
# @time W,H = nmf_euc(V,W)
# @time W,H = nmf_kl(V,W)
@time H = nmf_is(V,W)

# 結果を確認
println("Template NMF Result : $(norm(V - W*H))")

# NMF (Non-negative Matrix Factorization)
# Julia実装版

module NMF

export nmf_euc, nmf_kl, nmf_is

# EUC基準NMF処理（通常版）
# V     分解対象の行列
# k     因子分解数
# iter  イテレーション回数
function nmf_euc(V::AbstractMatrix,k::Int64=3;iter::Int64=100)
    # 変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative"
    @assert k > 0 "Factor Num : k must be greater than 0"
    @assert iter > 0 "Iteration Num : iter must be greater than 0"

    # W,Hを指定の形で乱数作成
    W = rand(size(V)[1],k)
    H = rand(k,size(V)[2])

    # 反復乗法更新
    for i in 1:iter
        W = W .* (V * H') ./ (W * H * H')
        H = H .* (W' * V) ./ (W' * W * H)
    end
    return W,H
end


# EUC基準NMF処理（テンプレートあり版）
# V     分解対象の行列
# W     テンプレート辞書行列
# iter  イテレーション回数
function nmf_euc(V::AbstractMatrix,W::AbstractMatrix,;iter::Int64=100)
    #変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative."
    @assert countnz(W.< 0) == 0 "Input Matrix W is not non-negative."
    @assert iter > 0 "Iteration Num : iter must be greater than 0."
    @assert size(V,1) == size(W,1) "Dimensions of V:$(size(V)) and W:$(size(W)) do not match. Align the first dimention."
    @assert size(W,2) > 0 "Factor Num : size(W,2) is greater than 0."

    # 指定の形でHを作成
    H = rand(size(W,2),size(V)[2])

    # H だけを更新
    for i in 1:iter
        H = H .* (W' * V) ./ (W' * W * H)
    end
    return H
end


# KL基準NMF処理（通常版）
# V     分解対象の行列
# k     因子分解数
# iter  イテレーション回数
function nmf_kl(V::AbstractMatrix,k::Int64=3;iter::Int64=100)
    #変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative"
    @assert k > 0 "Factor Num : k must be greater than 0"
    @assert iter > 0 "Iteration Num : iter must be greater than 0"

    # W,Hを指定の形で乱数作成
    W = rand(size(V)[1],k)
    H = rand(k,size(V)[2])

    # 反復乗法更新
    for i in 1:iter
        W = W .* ((V ./ (W*H)) * H') ./ sum(H',1)
        H = H .* (W' * (V ./ (W*H))) ./ sum(W',2)
    end
    return W,H
end


# KL基準NMF処理（テンプレートあり版）
# V     分解対象の行列
# W     テンプレート辞書行列
# iter  イテレーション回数
function nmf_kl(V::AbstractMatrix,W::AbstractMatrix;iter::Int64=100)
    #変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative."
    @assert countnz(W.< 0) == 0 "Input Matrix W is not non-negative."
    @assert iter > 0 "Iteration Num : iter must be greater than 0."
    @assert size(V,1) == size(W,1) "Dimensions of V:$(size(V)) and W:$(size(W)) do not match. Align the first dimention."
    @assert size(W,2) > 0 "Factor Num : size(W,2) is greater than 0."

    # 指定の形でHを作成
    H = rand(size(W,2),size(V)[2])

    # H だけを更新
    for i in 1:iter
        H = H .* (W' * (V ./ (W*H))) ./ sum(W',2)
    end
    return H
end


# IS基準NMF処理（通常版）
# V     分解対象の行列
# k     因子分解数
# iter  イテレーション回数
function nmf_is(V::AbstractMatrix,k::Int64=3;iter::Int64=100)
    #変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative"
    @assert k > 0 "Factor Num : k must be greater than 0"
    @assert iter > 0 "Iteration Num : iter must be greater than 0"

    # W,Hを指定の形で乱数作成
    W = rand(size(V)[1],k)
    H = rand(k,size(V)[2])

    # 反復乗法更新
    for i in 1:iter
        W = W .* sqrt.(((V ./ (W*H))*(H./sum(W*H,1))') ./(sum((H./sum(W*H,1))',1)))
        H = H.*sqrt.(((W ./ sum(W*H,2))' * (V./(W*H)) ./ (sum((W ./ sum(W*H,2))',2))))
    end
    return W,H
end


# IS基準NMF処理（テンプレートあり版）
# V     分解対象の行列
# W     テンプレート辞書行列
# iter  イテレーション回数
function nmf_is(V::AbstractMatrix,W::AbstractMatrix;iter::Int64=100)
    #変数チェック処理
    @assert countnz(V.< 0) == 0 "Input Matrix V is not non-negative."
    @assert countnz(W.< 0) == 0 "Input Matrix W is not non-negative."
    @assert iter > 0 "Iteration Num : iter must be greater than 0."
    @assert size(V,1) == size(W,1) "Dimensions of V:$(size(V)) and W:$(size(W)) do not match. Align the first dimention."
    @assert size(W,2) > 0 "Factor Num : size(W,2) is greater than 0."

    # 指定の形でHを作成
    H = rand(size(W,2),size(V)[2])

    # H だけを更新
    for i in 1:iter
        H = H.*sqrt.(((W ./ sum(W*H,2))' * (V./(W*H)) ./ (sum((W ./ sum(W*H,2))',2))))
    end
    return H
end

end

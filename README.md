# NMF_for_Julia
## Overview
JuliaでのNMF実装．

３種類の距離基準（EUC, KL, IS）のNMFを実装しています．

また，標準的な入力行列に対する行列分解を行う関数と，辞書行列にテンプレートを設定して励起行列のみを更新する形の関数も用意しています．

コード内の行列演算に問題がある可能性があるので，もしお気づきになられたら修正 or 連絡をお願い致します. 


## Usage
srcディレクトリの「Example.jl」に使い方を記載しています．


## Requirement
Julia 0.6.2


## License
Copyright © 2018 T_Sumida Distributed under the MIT License.

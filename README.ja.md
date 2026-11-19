# hetfamily — 異質な処置効果（jamovi モジュール）

ランダム化実験で「誰がより恩恵を受けるか」を調べます。汎用機械学習推測（Chernozhukov、Demirer、Duflo、Fernandez-Val）、優先度曲線、効果プロファイル、処置ルール。 詳しい英語の説明は [README.md](README.md)、中国語の要約は [README.zh.md](README.zh.md) を参照してください。

## 分析（jamovi メニュー：Heterogeneity）

| 分析 | 内容 |
|---|---|
| **異質な効果（BLP、GATES、CLAN）** | 最良線形予測器、予測効果で並べたグループ平均効果、最も影響を受けるグループと受けないグループの共変量平均 |
| **優先度曲線（RATE）** | TOC 曲線、AUTOC、Qini。ブートストラップ標準誤差 |
| **効果プロファイル（CATE 回帰）** | 二重頑健スコアを調整変数に回帰：平均での効果と標準偏差あたりの変化 |
| **処置ルール** | 予測効果がコストを上回る人だけを処置した場合の利得 |

各分析の出力には「方法の説明」と「データと設定ガイド」があります。

学習器：lasso、リッジ、弾性ネット、ランダムフォレスト、線形回帰。サンプル分割を多数回行い中央値で統合します。処置はランダム化されている（または傾向スコアが既知）必要があります。第 2 段階は `GenericML`・`grf` と一致を確認。サンプルデータ：`sim_het_trial.csv`、`sim_homogeneous_trial.csv`。

## インストール
- **jamovi ライブラリ**：jamovi 右上のモジュールボタンから *jamovi library* を開き、検索してインストール（モジュール公開後）。
- **サイドロード**：`.jmo` ファイルをダウンロードし、同じメニューの *Install from file* から導入します。`.jmo` はビルドした OS、プロセッサ、jamovi のバージョンでのみ動作します。
- **ソースから**：R と `jmvtools` が必要です：`jmvtools::install(home = "C:/Program Files/jamovi 28.17.0.0")`（パスはお使いの jamovi に合わせてください）。

## 引用とライセンス
引用：Sun, M. *hetfamily* (jamovi module). https://github.com/mingchao-sun/hetfamily 。使用した手法の原論文（出力に記載）も併せて引用してください。

ライセンス：`DESCRIPTION` を参照（GPL）。第三者のコードとデータは `NOTICE` に記載。不具合や提案：https://github.com/mingchao-sun/hetfamily/issues

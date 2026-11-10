# HGTFinder v0.1.0 — 公開用リリースノート案

GTDBタンパク質相同性検索のTop100から分類距離平均を計算し、最良ヒットが同属なら候補から除外する、探索的HGTスクリーニングCLIの初回公開予定版です。

- Python 3.11以上、検索はローカルDIAMOND。
- 同属最良ヒット除外を既定で有効化。保留判定・距離平均・採用ヒット・分母は維持。
- 検索済み結果の再集計、全query出力、ヒット監査、入力・コード・設定manifest。
- GTDB公式全種代表の分割参照DBに対応。DB本体は同梱しない。
- 合成fixtureとオフラインテスト、Python 3.11–3.13用CI定義。

LKSZ1事例では主要nif 5/5、全ゲノム1,148/3,285 CDSを候補化。精度ベンチマークやHGT確定を意味しない。同属間HGT、参照DB偏り、保存遺伝子、パラログに制約がある。

公開先：https://github.com/leoken422/HGTFinder 。MIT License、Copyright (c) 2026 Leonardo Ken Okumura。公開日は公開時に記録する。

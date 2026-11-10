# HGTFinder v0.1

GTDB由来の全CDSタンパク質参照DBを用い、単一原核生物ゲノムの全タンパク質からHGT候補を探索するPython CLIです。Python 3.11以上、検索には別途DIAMONDが必要です。Python実行時の外部依存はありません。

開発者：Leonardo Ken Okumura。MIT License（[LICENSE](LICENSE)）。公開先：https://github.com/leoken422/HGTFinder 。

**v0.1の標準条件：分類距離平均が閾値以上でも、採用済みTop100の最良ヒットが同属なら候補から除外します。** 探索的スクリーニングであり、HGT確率・統計的有意性・供与体確定を出しません。

## インストール

リポジトリルートで実行します。DB・入力配列は同梱していません。

```sh
python -m venv .venv
. .venv/bin/activate
python -m pip install .
hgtfinder --help
```

## 判定方法

1. 候補プールは最大1,000タンパク質/query。genome_idで自己を除外し、E値≤1e-5、query/subject coverage各50%以上を適用します。
2. 同一subjectを重複排除し、bit score降順、E値昇順、subject ID昇順で最大100タンパク質を選びます。完全同順位の重複行はalignment length・identity・coverage・座標等で決定的に並べます。
3. 選択100件内で同一ゲノムの最良1件を1票とし、計算可能な分類距離点の平均を主スコアにします。bit scoreの重みは掛けません。
4. 近縁参照充足性・最小10有効ヒット・最小5有効独立ゲノムの条件を満たさない場合は`pending`です。満たして平均≥2.0なら基本判定は`candidate`、未満なら`non_candidate`です。
5. 基本判定が候補で、選択済みTop100の最良1件が同属（同種を含む）なら`non_candidate`、理由`best_hit_same_genus`に変更します。保留や比較不能な属をこの条件で非候補には変更しません。

最良とは配列一致率最大でなくbit score最大です。同属ヒットが1件でもあれば除外するルールではありません。距離平均・ヒット集合・有効分母はこの条件で変わりません。

|対象との関係|点|
|---|---:|
|同一speciesの別ゲノム|0|
|同属・異種|1|
|同科・異属|2|
|同目・異科|3|
|同綱・異目|4|
|同門・異綱|5|
|同domain・異門|6|
|異domain|7|

自己は0点を記録して主集計から除外。距離を確定できないtaxonomy欠損はNAです。分類群は既知の祖先階級を含めて比較します。

## 検索せず計算を試す

同梱例は架空の配列・taxonomyです。ネットワーク、検索、DB構築を実行しません。

```sh
hgtfinder debug-score --hits examples/debug_hits.tsv \
  --metadata examples/reference_metadata.tsv \
  --query-taxonomy examples/query_taxonomy.tsv
```

q1は分子34、分母11、平均約3.091で基本判定は候補ですが、最良ヒットが同属なので最終判定は非候補です。期待出力は`examples/debug_expected.json`です。

## 入力と実行

- query FASTA：重複しないquery ID、全アミノ酸配列。
- query taxonomy TSV：`sample_id`, `genome_id`, `gtdb_release`, `gtdb_lineage`の単一サンプル行。genome_id不明なら空欄可ですが自己除外できません。
- reference FASTA：一意なsubject ID、各ゲノムの全CDSタンパク質。
- reference metadata TSV：`subject_id`, `genome_id`, `gtdb_lineage`。参照FASTAとID集合が完全一致する必要があります。
- 任意の遺伝子注釈TSV：`query_id`とgene/product/座標など。

GTDB lineageはd__/p__/c__/o__/f__/g__/s__の7階級。空階級は接頭辞のみ。NCBI分類やGTDB-TkマーカーDBで代用しません。DB・queryのGTDB releaseは一致させます。

```sh
hgtfinder prepare-db --reference reference.faa --metadata reference_metadata.tsv \
  --release R232 --selection species_representatives --protein-scope all-cds \
  --provenance '取得元、release、全CDS取得方法、チェックサム' --out outputs/db

hgtfinder run --query query.faa --query-taxonomy query_taxonomy.tsv \
  --db outputs/db --genes genes.tsv --out outputs/run
```

保存済み検索だけを再集計する場合：

```sh
hgtfinder summarize --query query.faa --query-taxonomy query_taxonomy.tsv \
  --db outputs/db --hits outputs/run/search/diamond.tsv \
  --search-manifest outputs/run/search/search_manifest.json \
  --genes genes.tsv --out outputs/reclassified
```

出力先は新規または空のディレクトリが必要です。検索manifestを渡すことで入力・検索条件・検索結果の対応を検証します。外部ヒットをmanifestなしで読む場合は由来未検証フラグを付けます。

設定変更は`hgtfinder --config config/custom.toml summarize ...`のようにサブコマンドの前に指定します。標準設定は`config/v0.1.toml`、大規模参照用設定は`config/full_gtdb_r232.toml`です。両方で同属除外を有効にしています。比較実験で無効化する場合はコピーした設定の`exclude_best_hit_same_genus = false`を明示し、標準v0.1と区別してください。

生DIAMOND出力の列順：

```text
qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore
```

## 出力

`genes.tsv`は0ヒットを含む全query、最終判定、理由、距離平均・分子分母、基本判定`base_distance_status`、`best_hit_same_genus`、`genus_veto_applied`を含みます。`hits_top100.tsv`、`hits_audit.tsv`、`taxonomy_counts.tsv`、`summary.json/md`、`run_manifest.json`を保存します。TSV欠損はNA、JSONはnullです。

## 検証

```sh
PYTHONPATH=src python -m unittest discover -s tests -v
```

標準テストはオフラインで、検索・DB構築の統合テスト2件は明示的な実行指定がない限りskipします。既存解析の検索をテストが再開することはありません。DIAMOND統合テストを別途実行する場合のみ、`HGTFINDER_TEST_DIAMOND`に絶対パスを設定します。統合テストは合成データ専用です。

## 適用範囲と制約

LKSZ1では標準条件で主要nif 5/5、対象領域31/34 CDS、全ゲノム1,148/3,285 CDSが候補となりました。これは単一事例であり、検出精度の検証ではありません。詳しくは[事例評価](docs/LKSZ1_CASE_STUDY.md)を参照してください。

参照DBの分類群構成、Top100打ち切り、保存遺伝子、パラログは未補正。同属除外は同属間HGTや同属に広がったHGTも落とし得ます。近縁ヒット不在は遺伝子欠如、非候補は垂直遺伝の証明ではありません。組成・遺伝子樹・パンゲノム・16S解析は本版の機能に含みません。NCBIリモートBLAST機能は公開パッケージに含めません。

分割参照DBについては[大規模参照の説明](docs/FULL_GTDB.md)、公開状況は[公開準備メモ](docs/PUBLICATION_STATUS.md)を参照してください。

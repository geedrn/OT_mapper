# OT mapper

## 概要
- 目的: CRISPR-Cas9のgRNA（SpCas9/NGG）に対するオフターゲット（OT）候補を取得し、遺伝子アノテーション（エクソン/イントロン）を付与し、Primer-BLAST用のプライマー設計リンクを生成します。
- 構成: シェルスクリプト2本（検出・マッピング）+ Rスクリプト1本（Primer-BLASTリンク生成）。

## ワークフロー
- OT候補検出（`scripts/OT_detector.sh`）
  - 対話入力で20ntのgRNA配列（PAMなし）を受け取り、自動で`NGG`を付与。
  - GGGenome（hg38）から候補配列CSVを取得（プラス/マイナス鎖、seed長は8または12ntを選択）。
  - 取得CSVを加工し、OT候補のBEDファイル `analysis/OT/OT_candidate.bed` を作成。
- 遺伝子アノテーション付与（`scripts/OT_mapper.sh`）
  - `bedtools intersect` により `OT_candidate.bed` と GENCODE V39（UCSC由来の加工BED）を突き合わせ、
    エクソン/イントロンの属性を付与した結果 `analysis/OT_mapper_results/OT_mapped.tsv` を出力。
- プライマーリンク生成（`scripts/primer_generate.R`）
  - `OT_mapped.tsv` を読み込み、各座標に対するPrimer-BLASTのURLを生成し、
    `analysis/OT_mapper_results/OT_with_primer.tsv` に出力。

## 使い方
1) 依存関係を準備
```bash
conda create -n OT -c bioconda r-optparse bedtools wget
conda activate OT
```

2) 実行
```bash
bash main.sh
```
- 対話で以下を入力します:
  - Spacer（20nt, PAMなし）例: `GCTGAAGCACTGCACGCCGT`
  - Seed長（8 または 12）

3) 出力物
- `analysis/OT/OT_candidate.bed`: OT候補のBED
- `analysis/OT_mapper_results/OT_mapped.tsv`: アノテーション付与済みTSV
- `analysis/OT_mapper_results/OT_with_primer.tsv`: Primer-BLAST URL付きTSV

## 入出力の詳細
- 入力（対話）
  - gRNA（20nt, PAMなし）。PAMはスクリプト内で`NGG`固定。
  - seed長（8または12）。
- アノテーション参照データ
  - `data/` 配下に GENCODE V39 を基にしたUCSC由来の加工BED
    - `UCSC_exons_modif_canonical.bed`（エクソン）
    - `UCSC_introns_modif_canonical.bed`（イントロン）
- 中間・最終出力
  - `analysis/OT/OT_candidate.bed`
  - `analysis/OT_mapper_results/OT_mapped.tsv`
  - `analysis/OT_mapper_results/OT_with_primer.tsv`

## 注意点・既知の制約
- ネットワーク必須: OT検出でGGGenomeからCSVを`wget`します（インターネット接続が必要）。
- PAM固定: 現状`NGG`のみ（SpCas9前提）。他PAMは未対応。
- bedtoolsの仕様: アノテーション非該当の候補は`OT_mapped.tsv`に現れません。`analysis/OT/OT_candidate.bed`も併せてご確認ください。
- タイムスタンプ: 現在のスクリプトは結果にタイムスタンプを付与していません。出力先は `analysis/` 配下です（従来記載の`data/`ではありません）。
- 列の整合性: `primer_generate.R` は入力TSVの1列目に染色体（例: `chr1`）を期待します。`scripts/OT_mapper.sh` の`cut`指定（`-f 5,1-3,8-9`）により列順が想定と異なる可能性があります。Primer生成でエラーが出る場合、`OT_mapped.tsv`の1–3列が`chr/start/end`になっているかをご確認ください。

## 参考: 実行例入力
- Spacer: `GCTGAAGCACTGCACGCCGT`
- Seed: `12`

## ライセンス/謝辞
- 本ツールはRyoおよびGabrielにより開発されました。

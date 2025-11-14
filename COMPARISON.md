# 実装比較: Bash vs R vs Shiny

## オプション/パラメータの比較

### Bash (OT_detector.sh)
- `-s, --spacer` (必須): 20nt gRNAスペーサー配列
- `-l, --seed-length` (デフォルト: 12): シード配列長 (8-12)
- `-p, --pam` (デフォルト: NGG): PAM配列
- `-g, --genome` (デフォルト: hg38): ゲノムアセンブリ
- `-f, --full-mismatch` (デフォルト: 3): 全長配列のミスマッチ許容数 (0-3)
- `-m, --seed-mismatch` (デフォルト: 1): シード配列のミスマッチ許容数 (0-3)
- `-o, --output-dir` (デフォルト: analysis): 出力ディレクトリ
- `--skip-annotation`: 遺伝子アノテーションをスキップ
- `--skip-primer`: Primer-BLASTリンク生成をスキップ

### R (analysis.R)
- `-s, --spacer` (必須): 20nt gRNAスペーサー配列
- `-l, --seed-length` (デフォルト: 12): シード配列長 (8-12)
- `-p, --pam` (デフォルト: NGG): PAM配列
- `-g, --genome` (デフォルト: hg38): ゲノムアセンブリ
- `-f, --full-mismatch` (デフォルト: 3): 全長配列のミスマッチ許容数 (0-3)
- `-m, --seed-mismatch` (デフォルト: 1): シード配列のミスマッチ許容数 (0-3)
- `-e, --exon-db`: エクソンアノテーションBEDファイル
- `-i, --intron-db`: イントロンアノテーションBEDファイル
- `-o, --output` (デフォルト: annotated_offtargets.tsv): 出力ファイル
- `--skip-primer`: Primer-BLASTリンク生成をスキップ

### Shiny (ui.R/server.R)
- `spacer` (テキスト入力): 20nt gRNAスペーサー配列
- `seed_length` (スライダー: 8-12, デフォルト: 12): シード配列長
- `full_mismatch` (スライダー: 0-3, デフォルト: 3): 全長配列のミスマッチ許容数
- `seed_mismatch` (スライダー: 0-3, デフォルト: 1): シード配列のミスマッチ許容数
- `pam` (選択: NGG, NRG, デフォルト: NGG): PAM配列
- `genome`: なし（固定でhg38）

## オプションの違い

### 一致しているオプション
✅ `spacer` (必須)
✅ `seed_length` (デフォルト: 12)
✅ `pam` (デフォルト: NGG)
✅ `full_mismatch` (デフォルト: 3)
✅ `seed_mismatch` (デフォルト: 1)

### 違いがあるオプション
❌ `genome`: BashとRでは選択可能、Shinyでは固定（hg38）
❌ `exon-db`, `intron-db`: Rのみ（Shinyではデフォルトパスを使用）
❌ `output-dir`/`output`: BashとRでは指定可能、Shinyでは自動
❌ `skip-annotation`: Bashのみ
❌ `skip-primer`: BashとRでは選択可能、Shinyでは常に実行

## アルゴリズムの比較

### 1. GGGenome API呼び出し
**すべて同一**
- URL形式: `https://gggenome.dbcls.jp/{genome}/{mismatch}/{strand}/nogap/{sequence}.csv`
- 全長配列とシード配列の両方を検索
- プラス鎖とマイナス鎖の両方を処理

### 2. PAMフィルタリング
**実装方法が異なる**

**Bash**:
- `sbjct`列（BEDファイルの4列目）から直接PAM部分を抽出
- プラス鎖: 末尾のPAM部分をチェック（例: `.GG$`）
- マイナス鎖: 先頭の逆相補PAM部分をチェック（例: `^CC.`）
- 正規表現パターンマッチングを使用

**R/Shiny**:
- `sbjct`列から直接PAM部分を抽出
- プラス鎖: `sbjct`の末尾からPAM部分を抽出
- マイナス鎖: `sbjct`の先頭からPAM部分を抽出し、逆相補パターンでマッチング
- 正規表現パターンマッチングを使用

**結論**: アルゴリズムは同一だが、実装の詳細が異なる

### 3. オーバーラップ検出
**すべて同一**
- `bedtools intersect`を使用
- コマンド: `bedtools intersect -b {seed_bed} -a {full_bed} -wa`
- 全長配列（`-a`）とシード配列（`-b`）の重なりを検出
- 全長配列を出力（`-wa`オプション）

**Bash**:
```bash
bedtools intersect -b "${temp_dir}/${strand}_seed_pam.bed" \
                   -a "${temp_dir}/${strand}_full_pam.bed" \
                   -wa > "${temp_dir}/${strand}_overlaps.bed"
```

**R**:
```r
cmd <- paste("bedtools intersect -b", plus_seed_bed, "-a", plus_full_bed, "-wa >", plus_overlap_bed)
system(cmd)
```

**結論**: 完全に同一

### 4. データ処理フロー
**すべて同一**
1. GGGenome APIからデータをダウンロード
2. PAMフィルタリング
3. オーバーラップ検出（bedtools intersect）
4. 結果の結合（プラス鎖 + マイナス鎖）
5. 遺伝子アノテーション（オプション）
6. Primer-BLASTリンク生成（オプション）

## まとめ

### アルゴリズム
✅ **同一**: すべての実装で同じアルゴリズムを使用
- GGGenome API呼び出し方法
- オーバーラップ検出（bedtools intersect）
- データ処理フロー

⚠️ **実装の詳細が異なる**: PAMフィルタリングの実装方法は異なるが、ロジックは同一

### オプション
✅ **主要オプションは一致**: spacer, seed_length, pam, full_mismatch, seed_mismatch
❌ **一部のオプションが異なる**: genome選択、出力先指定、スキップオプション

### 推奨事項
1. **Shinyに`genome`選択を追加**: BashとRとの一貫性のため
2. **PAMフィルタリングの実装を統一**: より確実な結果のため
3. **Shinyに`skip-primer`オプションを追加**: パフォーマンス向上のため


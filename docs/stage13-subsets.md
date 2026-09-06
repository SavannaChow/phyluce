# Stage 13：獨立樣本分支分析

完整流程 `f` 仍執行 stage 1–12。Stage 13 必須另外選擇。
每個新分支以原始專案的 stage 5 loci 紀錄與 stage 6 已擷取序列為來源，
依保留名單重新執行 stage 7–12，包括 MAFFT 對齊、修剪與完整度篩選。
原始 stage 1–12、`project.info` 和原本的 log 不會被修改。

## 互動操作

1. 啟動 `phyluce_launcher`，選擇 `13`，輸入原始專案的位置或名稱。
2. 選 `create`，輸入你要的分析名稱。
3. 打開新資料夾外面的 **`original_samples.txt`**，刪掉不要分析的樣本整行，再存檔。
   這份檔案雖然名為 original samples，**編輯後留下的就是本次要分析的樣本**。
   一行一個完整 sample ID，不要改寫 ID，也不要加入標題或物種名稱。
4. 再次選 stage `13` → `run`，依編號選擇該資料夾。
5. 第一次執行可調整 stage 7–12 參數，預設沿用建立分支時的原專案設定。
6. 資料準備完成後，執行該分支 stage 12 裡產生的 IQ-TREE／ASTRAL 腳本。

與原本 stage 12 一樣，stage 13 **產生建樹資料與腳本，不會自動啟動建樹**。
`prepared` 表示 stage 7–12 的資料準備完成，不表示 IQ-TREE 或 ASTRAL 已跑完。

產生的 IQ-TREE 腳本旁會附上 `watch_iqtree_progress_htop.sh`，執行分析時會
自動嘗試開啟另一個進度視窗。詳見 [IQ-TREE 進度監看](iqtree-progress.md)。

## 資料夾

日期時間使用台灣時間，格式為 `名稱_YYYYMMDD_HHMMSS`。
同一秒建立同名分支時會加 `_02` 等編號，避免覆寫。

```text
原始專案/
├── 01_… 到 12_Analysis_Branches/
└── 13_Subset_Analyses/
    ├── analysis_index.tsv
    └── 你的名稱_20260904_143025/
        ├── original_samples.txt          ← 只需編輯這份
        ├── analysis_records/
        │   ├── all_samples.tsv           ← 原始可用樣本與擷取 loci 數
        │   ├── analysis_samples.txt      ← 開始執行時保存的分析名單
        │   ├── excluded_samples_in_this_run.txt
        │   ├── analysis_info.json        ← 來源、參數、狀態及續跑紀錄
        │   ├── final_sample_presence.tsv ← 各最終矩陣實際包含的樣本
        │   └── run.lock
        ├── project.info                  ← 分支專用的輸出識別名稱
        ├── logs/
        ├── 07_Aggregated_loci/
        ├── 08_Export_Locus_FASTAs/
        ├── 09_Aligned_Locus_FASTAs/
        ├── 10_Trimmed_Locus_ALIGNMENTS/
        ├── 11_Filtered_Locus_ALIGNMENTS/
        └── 12_Analysis_Branches/
```

執行紀錄、log 和 stage 目錄會在需要時產生。各分支不複製 stage 1–6。
`all_samples.tsv` 根據原始 stage 6 可用序列建立，不是從最後一棵樹取得名單。
`analysis_samples.txt` 記錄本次選入的樣本；若部分樣本經 loci 篩選後沒有進入最終矩陣，
程式會顯示訊息，並記錄在 `final_sample_presence.tsv`。

## 多次分析與續跑

- `create`：每次都先列出原始專案的全部可用樣本，各次排除不會累積。
- `copy`：選擇既有分支，再輸入新名稱。只複製樣本選擇與參數，仍直接使用原始資料。
  已開始的分支以保存的分析名單為準；草稿以可編輯名單為準。
  想加回樣本，可從 `analysis_records/all_samples.tsv` 複製完整 ID 到新名單。
- `run`：草稿開始分析；失敗或中斷的分支從最後完成的 stage 之後繼續。
  尚未完成的 stage 會清除其分支內的部分輸出後重跑。
- 已完成的分支再次執行 `run` 時，會檢查已記錄完成的 Stage 7–12 資料夾。
  如果某個 Stage 資料夾已被整個刪除或是空的，會從該 Stage 開始重跑，並重建其後續 Stage；
  例如只刪除 `12_Analysis_Branches/`，就只重建 Stage 12，保留 Stage 7–11。
- `list`：列出分支和狀態。`analysis_index.tsv` 另外提供樣本數、完整度門檻及進度。

開始後不能更改同一分支的樣本名單或參數再續跑；應 `copy` 成新分支。
程式會拒絕重複／不存在的 ID、空名單，以及少於三個樣本的組合。
所選樣本的 stage 5/6 來源會比對建立分支時的 SHA-256；來源已改變時，
請建立新分支，以免把不同批次的資料混進已開始的分析。

完整度使用新的樣本數計算，整數取整規則沿用原 stage 11。
本功能沿用 stage 5/6 的位點判定與擷取方式，不重新搜尋 UCE 或改動原始序列。

## 指令模式

以下從專案上一層執行，將 `/path/to/phyluce_launcher` 換成腳本實際位置。

```bash
# 建立名單；此時尚不執行對齊
python3 /path/to/phyluce_launcher --project-id MyProject --start-stage 13 \
  --subset-action create --subset-name ingroup_only

# 編輯 original_samples.txt 後，使用實際建立的資料夾名稱
python3 /path/to/phyluce_launcher --project-id MyProject --start-stage 13 \
  --subset-action run --subset-dir ingroup_only_20260904_143025

# 只預覽，不寫入分析結果或執行紀錄
python3 /path/to/phyluce_launcher --project-id MyProject --start-stage 13 \
  --subset-action run --subset-dir ingroup_only_20260904_143025 --dry-run

# 複製成另一個獨立組合
python3 /path/to/phyluce_launcher --project-id MyProject --start-stage 13 \
  --subset-action copy --subset-dir ingroup_only_20260904_143025 \
  --subset-name ingroup_with_outgroup
```

首次 `run` 可使用既有的 `--dataset-modes`、`--post-align-thresholds`、
`--trim-mode`、`--mafft-threads-per-job` 等參數；續跑會使用已保存的設定。

## 開發驗證

`python3 -m pytest phyluce/tests/test_subsets.py` 檢查分支隔離、名單、
重新篩選完整度、來源檔案不變、複製及中斷續跑。
需 pytest、Biopython 與 NumPy；測試以預先對齊的小型序列搭配 MAFFT 替身，
不驗證真實 MAFFT 的對齊品質，也不執行 IQ-TREE／ASTRAL。

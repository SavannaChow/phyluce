# IQ-TREE 分析與進度視窗

Stage 12 產生的 supermatrix 與 gene-tree IQ-TREE 分析腳本，會在同一個資料夾附上
`watch_iqtree_progress_htop.sh`。Stage 13 各分支使用同一套 stage 12 輸出方式，因此也適用。
這是終端機進度畫面，不需要另外安裝 htop。

直接執行產生的 `run_iqtree_*.sh` 或 `01_run_iqtree_gene_trees.sh`：

- IQ-TREE 在目前終端執行，進度監看另外開啟。
- Supermatrix 自動使用該次 `--prefix` 對應的 `.log`，不用手動找檔案。
- Gene-tree 只開一個總進度畫面，顯示完成的 loci 數量與錯誤 log 數量。
- Log 尚未產生時會等待；分析結束或失敗時顯示最後狀態並結束監看。
- 關閉監看視窗或按 Ctrl-C 只停止監看，不會停止 IQ-TREE。

## 開啟方式

若已在 tmux／screen 內，開啟新的視窗；macOS 本機使用 Terminal；
Linux 桌面使用可用的 gnome-terminal、konsole、x-terminal-emulator 或 xterm。
Linux 圖形終端啟動後仍會印出手動指令，以便桌面環境拒絕開窗時使用。

在無桌面的 SSH 環境，若有 tmux，會建立獨立的監看 session，並印出
另一個 SSH 終端可以執行的 `tmux attach -t ...` 指令。
若沒有可用的視窗工具，分析照常執行，終端會印出完整的手動監看指令。
無法透過遠端 shell 直接開啟你本機的終端視窗。

## 選用設定

```bash
# 每 10 秒更新進度（預設 60 秒）
IQTREE_PROGRESS_INTERVAL=10 bash run_iqtree_你的分析名稱.sh

# 本次不自動開啟監看
IQTREE_PROGRESS=0 bash run_iqtree_你的分析名稱.sh

# 隨時在另一個終端手動監看指定 log
bash watch_iqtree_progress_htop.sh 你的分析.log 10
```

分析腳本另外寫入 `.progress_status` 狀態檔：supermatrix 使用
`<prefix>.progress_status`，gene-tree 使用 `iqtree_gene_trees/.progress_status`。
自動啟動的監看會讀取它，以區分正常完成和分析失敗。
Gene-tree 中原本已存在 `.treefile` 而被跳過的 loci，依現有分析腳本的規則處理；
整批成功結束時，總進度顯示全部完成。

此功能適用於更新後重新產生的腳本。已經輸出的舊分析資料夾不會自動更動；
須重新執行該組合的 stage 12 資料準備，才會得到新版分析腳本及監看腳本。
已在執行的 IQ-TREE 工作不會因更新程式而被中斷。

ASTRAL 腳本本身不啟動 IQ-TREE 監看；gene-tree 分支的監看是在
`01_run_iqtree_gene_trees.sh` 執行時開啟。

#!/usr/bin/env bash
set -u

LOG_FILE="${1:-}"
INTERVAL="${2:-60}"

if [[ -z "$LOG_FILE" ]]; then
    LOG_FILE="$(find . -maxdepth 1 -type f -name '*.log' -print -quit)"
fi

if [[ -z "$LOG_FILE" || ! -f "$LOG_FILE" ]]; then
    echo "用法：$0 IQTREE.log [更新秒數]"
    echo "例如：$0 AhyaTW_edge-incomplete-min_taxa_050.log 60"
    exit 1
fi

format_seconds() {
    local s="${1:-0}"
    (( s < 0 )) && s=0
    printf '%02dh %02dm %02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
}

# 回傳符合 pattern 的最後一個行號；找不到時回傳 0
last_line_no() {
    local pattern="$1"
    local n
    n="$(grep -nE "$pattern" "$LOG_FILE" 2>/dev/null | tail -1 | cut -d: -f1 || true)"
    echo "${n:-0}"
}

while true; do
    # 各階段最後出現的位置
    line_model="$(last_line_no 'In ModelFinder|ModelFinder requires|Model selection')"
    line_merge="$(last_line_no 'Merging partitions|Partition merging|merging partitions|merge partitions')"
    line_model_done="$(last_line_no 'CPU time for ModelFinder|Wall-clock time for ModelFinder')"
    line_ufboot="$(last_line_no 'Generating [0-9]+ samples for ultrafast bootstrap|Ultrafast bootstrap|Starting.*bootstrap')"
    line_bootstrap_support="$(last_line_no 'Bootstrap support|Bootstrap correlation coefficient')"
    line_finished="$(last_line_no 'Analysis finished|Date and Time:.*finished|Total wall-clock time used')"

    # 預設階段
    phase="初始化／讀取資料"

    # 依「最後出現位置」判斷目前階段。
    # 關鍵：一旦出現 Generating N samples for ultrafast bootstrap，
    # 後續 Iteration 就屬於 UFBoot，不再被前面的 MERGE 誤判。
    if (( line_model > 0 )); then
        phase="ModelFinder／模型選擇"
    fi

    if (( line_merge > line_model )) && (( line_merge > line_ufboot )); then
        phase="partition model 合併（MERGE）"
    fi

    if (( line_model_done > line_merge )) && (( line_ufboot == 0 )); then
        phase="ModelFinder／MERGE 已完成，準備樹搜尋／bootstrap"
    fi

    if (( line_ufboot > 0 )) && (( line_ufboot > line_merge )); then
        phase="Ultrafast Bootstrap"
    fi

    if (( line_bootstrap_support > line_ufboot )) && (( line_ufboot > 0 )); then
        phase="Ultrafast Bootstrap"
    fi

    if (( line_finished > line_ufboot )) && (( line_finished > 0 )); then
        phase="輸出結果／收尾"
    fi

    # 找最近一筆 iteration
    latest="$(grep -E 'Iteration [0-9]+ / LogL:' "$LOG_FILE" 2>/dev/null | tail -1 || true)"

    # 狀態
    status="執行中"
    if (( line_finished > 0 )); then
        status="已完成"
    fi

    # 最近關鍵訊息
    key_line="$(
        grep -E \
        'In ModelFinder|CPU time for ModelFinder|Wall-clock time for ModelFinder|Merging partitions|Partition merging|merging partitions|merge partitions|Generating [0-9]+ samples for ultrafast bootstrap|Bootstrap correlation coefficient|Iteration [0-9]+ / LogL:|Analysis finished|Total wall-clock time used|BEST SCORE' \
        "$LOG_FILE" 2>/dev/null | tail -1 || true
    )"

    clear
    echo "IQ-TREE 進度監控：$(date '+%F %T')"
    echo "Log：$LOG_FILE"
    echo "狀態：$status"
    echo "目前階段：$phase"

    # 只有 UFBoot 才把 Iteration 當 bootstrap replicate
    if [[ "$phase" == "Ultrafast Bootstrap" ]] && [[ -n "$latest" ]]; then
        iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$latest")"

        # 從 log 自動讀取 UFBoot replicate 總數，不硬寫 1000
        bootstrap_total="$(
            grep -E 'Generating [0-9]+ samples for ultrafast bootstrap' "$LOG_FILE" 2>/dev/null \
            | tail -1 \
            | sed -E 's/.*Generating ([0-9]+) samples for ultrafast bootstrap.*/\1/' || true
        )"

        elapsed="$(sed -E 's/.*Time: ([0-9]+)h:([0-9]+)m:([0-9]+)s.*/\1 \2 \3/' <<< "$latest")"
        read -r eh em es <<< "$elapsed"
        elapsed_sec=$((10#$eh*3600 + 10#$em*60 + 10#$es))

        eta="$(sed -E 's/.*\(([0-9]+)h:([0-9]+)m:([0-9]+)s left\).*/\1 \2 \3/' <<< "$latest")"
        if [[ "$eta" =~ ^[0-9]+[[:space:]][0-9]+[[:space:]][0-9]+$ ]]; then
            read -r rh rm rs <<< "$eta"
            iqtree_remaining=$((10#$rh*3600 + 10#$rm*60 + 10#$rs))
        else
            iqtree_remaining=-1
        fi

        if [[ "$bootstrap_total" =~ ^[0-9]+$ ]] && (( bootstrap_total > 0 )); then
            percent=$((iteration * 100 / bootstrap_total))

            if (( iteration > 0 && iteration < bootstrap_total )); then
                avg_remaining=$(( elapsed_sec * (bootstrap_total - iteration) / iteration ))
            else
                avg_remaining=0
            fi

            echo "Bootstrap 進度：$iteration / $bootstrap_total (${percent}%)"
            echo "Bootstrap 已用時間：$(format_seconds "$elapsed_sec")"
            echo "依目前總平均速度估算剩餘：$(format_seconds "$avg_remaining")"

            if (( iqtree_remaining >= 0 )); then
                echo "IQ-TREE 最近顯示剩餘：$(format_seconds "$iqtree_remaining")"
            fi
        else
            echo "Bootstrap iteration：$iteration"
            echo "Bootstrap 已用時間：$(format_seconds "$elapsed_sec")"
            if (( iqtree_remaining >= 0 )); then
                echo "IQ-TREE 最近顯示剩餘：$(format_seconds "$iqtree_remaining")"
            fi
        fi

    elif [[ -n "$latest" ]]; then
        # 非 UFBoot 階段只顯示 iteration，不假裝它是 /1000
        iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$latest")"
        echo "目前 Iteration：$iteration"
    else
        echo "目前尚未找到 Iteration 記錄"
    fi

    echo
    if [[ -n "$latest" ]]; then
        echo "最近一筆：$latest"
    fi
    if [[ -n "$key_line" ]]; then
        echo "最近關鍵訊息：$key_line"
    fi
    echo
    echo "Ctrl-C 離開監控；不會停止 IQ-TREE。"

    [[ "$status" == "已完成" ]] && exit 0
    sleep "$INTERVAL"
done

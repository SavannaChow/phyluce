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

while true; do
    latest="$(grep -E 'Iteration [0-9]+ / LogL:' "$LOG_FILE" | tail -1)"
    if [[ -z "$latest" ]]; then
        clear
        echo "$(date '+%F %T')"
        echo "尚未找到 IQ-TREE iteration：$LOG_FILE"
        sleep "$INTERVAL"
        continue
    fi

    iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$latest")"
    elapsed="$(sed -E 's/.*Time: ([0-9]+)h:([0-9]+)m:([0-9]+)s.*/\1 \2 \3/' <<< "$latest")"
    read -r eh em es <<< "$elapsed"
    elapsed_sec=$((eh*3600 + em*60 + es))
    percent=$((iteration * 100 / 1000))
    remaining=$(( (elapsed_sec * (1000 - iteration)) / iteration ))
    eta="$(sed -E 's/.*\(([0-9]+)h:([0-9]+)m:([0-9]+)s left\).*/\1 \2 \3/' <<< "$latest")"
    read -r rh rm rs <<< "$eta"
    iqtree_remaining=$((rh*3600 + rm*60 + rs))

    if grep -q 'Iteration 1000 /' "$LOG_FILE"; then
        status="已完成 1000 / 1000"
    else
        status="執行中"
    fi

    phase="樹搜尋／bootstrap"
    if grep -qE 'ModelFinder|Model selection' "$LOG_FILE"; then
        phase="ModelFinder／模型選擇"
    fi
    if grep -qE 'Merging partitions|Partition merging|MERGE' "$LOG_FILE"; then
        phase="partition model 合併（MERGE）"
    fi
    if grep -qE 'Bootstrap analysis|Ultrafast bootstrap|Starting.*bootstrap|Bootstrap support' "$LOG_FILE"; then
        phase="Ultrafast bootstrap"
    fi
    if grep -qE 'Analysis finished|IQ-TREE multicore version|Date and Time' "$LOG_FILE" && [[ "$iteration" -ge 1000 ]]; then
        phase="輸出結果／收尾"
    fi
    key_line="$(grep -E 'ModelFinder|Model selection|Merging partitions|Partition merging|Bootstrap|Iteration|Analysis finished|BEST SCORE' "$LOG_FILE" | tail -1)"

    clear
    echo "IQ-TREE 進度監控：$(date '+%F %T')"
    echo "Log：$LOG_FILE"
    echo "狀態：$status"
    echo "目前階段：$phase"
    echo "進度：$iteration / 1000 (${percent}%)"
    echo "已用時間：$(format_seconds "$elapsed_sec")"
    echo "依總平均速度估算剩餘：$(format_seconds "$remaining")"
    echo "IQ-TREE 最近顯示剩餘：$(format_seconds "$iqtree_remaining")"
    echo
    echo "最近一筆：$latest"
    echo "最近關鍵訊息：$key_line"
    echo
    echo "Ctrl-C 離開監控；不會停止 IQ-TREE。"

    [[ "$iteration" -ge 1000 ]] && exit 0
    sleep "$INTERVAL"
done

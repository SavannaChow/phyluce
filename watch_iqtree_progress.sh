#!/usr/bin/env bash
set -u

LOG_FILE="${1:-}"
INTERVAL="${2:-60}"
if [[ -z "$LOG_FILE" || ! -f "$LOG_FILE" ]]; then
    echo "用法：$0 IQTREE.log [更新秒數]"; exit 1
fi

format_seconds() {
    local s="${1:-0}"; (( s < 0 )) && s=0
    printf '%02dh %02dm %02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
}
file_size() { wc -c < "$LOG_FILE" | tr -d ' '; }

# 僅第一次讀取既有 log 尾端；之後按 byte offset 增量讀取，不會漏訊息。
offset=0; phase="初始化／讀取資料"; last_iteration=""
last_iteration_line=""; bootstrap_total=""; last_key_line=""; finished=0; pending=""

while true; do
    size="$(file_size)"
    if (( size < offset )); then offset=0; phase="初始化／讀取資料"; fi
    if (( offset == 0 && size > 200000 )); then
        chunk="$(tail -c 200000 "$LOG_FILE")"; offset="$size"
    elif (( size > offset )); then
        chunk="$(tail -c +$((offset + 1)) "$LOG_FILE")"; offset="$size"
    else
        chunk=""
    fi

    # 保留尚未寫完的最後一行，避免跨兩次輪詢時漏掉關鍵字。
    if [[ -n "$chunk" ]]; then
        data="${pending}${chunk}"
        if [[ "$data" == *$'\n'* ]]; then
            pending="${data##*$'\n'}"
            chunk="${data%$'\n'*}"
        else
            pending="$data"
            chunk=""
        fi
    fi

    if [[ -n "$chunk" ]]; then
        grep -qE 'ModelFinder|Model selection' <<< "$chunk" && phase="ModelFinder／模型選擇"
        grep -qEi 'merg(e|ing) partitions|partition merging' <<< "$chunk" && phase="partition model 合併（MERGE）"
        grep -qE 'STARTING TREE SEARCH|Starting tree search|Performing tree search|BEST SCORE|Optimizing NNI' <<< "$chunk" && phase="ML 樹搜尋／參數最佳化"
        grep -qEi 'SH-aLRT|aLRT|Approximate likelihood-ratio' <<< "$chunk" && phase="SH-aLRT／分支支持度計算"
        grep -qEi 'ultrafast bootstrap|Generating [0-9]+ samples' <<< "$chunk" && phase="Ultrafast Bootstrap"
        grep -qE 'Analysis finished|Total wall-clock time used' <<< "$chunk" && { phase="輸出結果／收尾"; finished=1; }

        line="$(grep -E 'Iteration [0-9]+ / LogL:' <<< "$chunk" | tail -1 || true)"
        if [[ -n "$line" ]]; then
            last_iteration_line="$line"
            last_iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$line")"
        fi
        found_total="$(grep -Ei 'Generating [0-9]+ samples for ultrafast bootstrap' <<< "$chunk" | tail -1 | sed -E 's/.*Generating ([0-9]+) samples.*/\1/' || true)"
        [[ "$found_total" =~ ^[0-9]+$ ]] && bootstrap_total="$found_total"
        key="$(grep -E 'ModelFinder|Merging partitions|Partition merging|STARTING TREE SEARCH|Starting tree search|SH-aLRT|aLRT|ultrafast bootstrap|Iteration [0-9]+ / LogL:|Analysis finished|Total wall-clock time used' <<< "$chunk" | tail -1 || true)"
        [[ -n "$key" ]] && last_key_line="$key"
    fi

    clear
    echo "IQ-TREE 進度監控：$(date '+%F %T')"
    echo "目前階段：$phase"
    echo "狀態：$([[ $finished == 1 ]] && echo 已完成 || echo 執行中)"
    if [[ "$phase" == "Ultrafast Bootstrap" && "$last_iteration" =~ ^[0-9]+$ ]]; then
        total="${bootstrap_total:-1000}"
        echo "Bootstrap 進度：$last_iteration / $total ($((last_iteration * 100 / total))%)"
        elapsed="$(sed -nE 's/.*Time: ([0-9]+)h:([0-9]+)m:([0-9]+)s.*/\1 \2 \3/p' <<< "$last_iteration_line")"
        if [[ -n "$elapsed" ]]; then
            read -r h m s <<< "$elapsed"; elapsed_sec=$((10#$h*3600 + 10#$m*60 + 10#$s))
            remain=$(( elapsed_sec * (total-last_iteration) / last_iteration ))
            echo "Bootstrap 預估剩餘：$(format_seconds "$remain")"
        fi
    elif [[ -n "$last_iteration" ]]; then echo "目前 Iteration：$last_iteration"; fi
    [[ -n "$last_key_line" ]] && echo "最近關鍵訊息：$last_key_line"
    echo "每 ${INTERVAL} 秒增量檢查一次；Ctrl-C 不會停止 IQ-TREE。"
    [[ $finished == 1 ]] && exit 0
    sleep "$INTERVAL"
done

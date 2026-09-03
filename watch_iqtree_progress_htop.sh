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

# ANSI 顏色
RESET=$'\033[0m'
BOLD=$'\033[1m'
DIM=$'\033[2m'
RED=$'\033[31m'
GREEN=$'\033[32m'
YELLOW=$'\033[33m'
BLUE=$'\033[34m'
MAGENTA=$'\033[35m'
CYAN=$'\033[36m'
WHITE=$'\033[37m'

# 隱藏游標，離開時恢復
printf '\033[?25l'
trap 'printf "\033[?25h\033[0m\n"' EXIT INT TERM

format_seconds() {
    local s="${1:-0}"
    (( s < 0 )) && s=0
    printf '%02dh %02dm %02ds' $((s/3600)) $(((s%3600)/60)) $((s%60))
}

last_line_no() {
    local pattern="$1"
    local n
    n="$(grep -nE "$pattern" "$LOG_FILE" 2>/dev/null | tail -1 | cut -d: -f1 || true)"
    echo "${n:-0}"
}

term_width() {
    local w
    w="$(tput cols 2>/dev/null || echo 100)"
    (( w < 80 )) && w=80
    (( w > 140 )) && w=140
    echo "$w"
}

repeat_char() {
    local char="$1"
    local count="$2"
    local i
    for ((i=0; i<count; i++)); do
        printf '%s' "$char"
    done
}

progress_bar() {
    local current="$1"
    local total="$2"
    local width="$3"

    if (( total <= 0 )); then
        printf "[%s]" "$(repeat_char " " "$width")"
        return
    fi

    (( current < 0 )) && current=0
    (( current > total )) && current=$total

    local filled=$(( current * width / total ))
    local empty=$(( width - filled ))

    printf "${GREEN}["
    (( filled > 0 )) && printf "%s" "$(repeat_char "█" "$filled")"
    (( empty > 0 )) && printf "%s" "$(repeat_char "░" "$empty")"
    printf "]${RESET}"
}

phase_color() {
    case "$1" in
        *"Ultrafast Bootstrap"*) echo "$GREEN" ;;
        *"MERGE"*) echo "$MAGENTA" ;;
        *"ModelFinder"*) echo "$CYAN" ;;
        *"收尾"*) echo "$YELLOW" ;;
        *"初始化"*) echo "$BLUE" ;;
        *) echo "$WHITE" ;;
    esac
}

while true; do
    line_model="$(last_line_no 'In ModelFinder|ModelFinder requires|Model selection')"
    line_merge="$(last_line_no 'Merging partitions|Partition merging|merging partitions|merge partitions')"
    line_model_done="$(last_line_no 'CPU time for ModelFinder|Wall-clock time for ModelFinder')"
    line_ufboot="$(last_line_no 'Generating [0-9]+ samples for ultrafast bootstrap|Ultrafast bootstrap|Starting.*bootstrap')"
    line_bootstrap_support="$(last_line_no 'Bootstrap support|Bootstrap correlation coefficient')"
    line_finished="$(last_line_no 'Analysis finished|Date and Time:.*finished|Total wall-clock time used')"

    phase="初始化／讀取資料"

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

    latest="$(grep -E 'Iteration [0-9]+ / LogL:' "$LOG_FILE" 2>/dev/null | tail -1 || true)"

    status="RUNNING"
    status_color="$GREEN"
    if (( line_finished > 0 )); then
        status="FINISHED"
        status_color="$CYAN"
    fi

    key_line="$(
        grep -E \
        'In ModelFinder|CPU time for ModelFinder|Wall-clock time for ModelFinder|Merging partitions|Partition merging|merging partitions|merge partitions|Generating [0-9]+ samples for ultrafast bootstrap|Bootstrap correlation coefficient|Iteration [0-9]+ / LogL:|Analysis finished|Total wall-clock time used|BEST SCORE' \
        "$LOG_FILE" 2>/dev/null | tail -1 || true
    )"

    W="$(term_width)"
    INNER=$((W - 4))
    BAR_WIDTH=$((W - 34))
    (( BAR_WIDTH < 20 )) && BAR_WIDTH=20
    (( BAR_WIDTH > 70 )) && BAR_WIDTH=70

    clear

    # Header
    printf "${BOLD}${CYAN}┌"
    repeat_char "─" $((W - 2))
    printf "┐${RESET}\n"

    title=" IQ-TREE Progress Monitor "
    printf "${BOLD}${CYAN}│${RESET}${BOLD} %-*s ${CYAN}│${RESET}\n" $((W - 4)) "$title"

    printf "${BOLD}${CYAN}├"
    repeat_char "─" $((W - 2))
    printf "┤${RESET}\n"

    printf "${CYAN}│${RESET} ${DIM}Time${RESET}   %-19s   ${DIM}Status${RESET}  ${status_color}${BOLD}%-10s${RESET}%*s${CYAN}│${RESET}\n" \
        "$(date '+%F %T')" "$status" $((W - 52)) ""

    phase_c="$(phase_color "$phase")"
    printf "${CYAN}│${RESET} ${DIM}Phase${RESET}  ${phase_c}${BOLD}%-*s${RESET} ${CYAN}│${RESET}\n" \
        $((W - 12)) "$phase"

    printf "${CYAN}│${RESET} ${DIM}Log${RESET}    %-*s ${CYAN}│${RESET}\n" \
        $((W - 12)) "$(basename "$LOG_FILE")"

    printf "${BOLD}${CYAN}├"
    repeat_char "─" $((W - 2))
    printf "┤${RESET}\n"

    if [[ "$phase" == "Ultrafast Bootstrap" ]] && [[ -n "$latest" ]]; then
        iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$latest")"

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

            printf "${CYAN}│${RESET} ${BOLD}Bootstrap progress${RESET}%*s${CYAN}│${RESET}\n" $((W - 22)) ""
            printf "${CYAN}│${RESET}   "
            progress_bar "$iteration" "$bootstrap_total" "$BAR_WIDTH"
            printf "  ${BOLD}%3d%%${RESET}  %d / %d%*s${CYAN}│${RESET}\n" \
                "$percent" "$iteration" "$bootstrap_total" \
                $((W - BAR_WIDTH - 26 - ${#iteration} - ${#bootstrap_total})) ""

            printf "${CYAN}│${RESET}   ${DIM}Elapsed${RESET}       ${BOLD}%-14s${RESET} ${DIM}Avg ETA${RESET}  ${BOLD}%-14s${RESET}%*s${CYAN}│${RESET}\n" \
                "$(format_seconds "$elapsed_sec")" "$(format_seconds "$avg_remaining")" $((W - 59)) ""

            if (( iqtree_remaining >= 0 )); then
                printf "${CYAN}│${RESET}   ${DIM}IQ-TREE ETA${RESET}   ${YELLOW}${BOLD}%-14s${RESET}%*s${CYAN}│${RESET}\n" \
                    "$(format_seconds "$iqtree_remaining")" $((W - 34)) ""
            fi
        else
            printf "${CYAN}│${RESET} ${BOLD}Bootstrap iteration${RESET}  %d%*s${CYAN}│${RESET}\n" \
                "$iteration" $((W - 27 - ${#iteration})) ""
        fi

    elif [[ -n "$latest" ]]; then
        iteration="$(sed -E 's/.*Iteration ([0-9]+).*/\1/' <<< "$latest")"
        printf "${CYAN}│${RESET} ${BOLD}Current iteration${RESET}  ${YELLOW}%d${RESET}%*s${CYAN}│${RESET}\n" \
            "$iteration" $((W - 25 - ${#iteration})) ""
        printf "${CYAN}│${RESET} ${DIM}This stage has no fixed total, so no percentage is shown.${RESET}%*s${CYAN}│${RESET}\n" \
            $((W - 62)) ""
    else
        printf "${CYAN}│${RESET} ${DIM}No Iteration record found yet.${RESET}%*s${CYAN}│${RESET}\n" \
            $((W - 34)) ""
    fi

    printf "${BOLD}${CYAN}├"
    repeat_char "─" $((W - 2))
    printf "┤${RESET}\n"

    printf "${CYAN}│${RESET} ${BOLD}Recent IQ-TREE message${RESET}%*s${CYAN}│${RESET}\n" $((W - 26)) ""

    if [[ -n "$latest" ]]; then
        short_latest="${latest:0:$((W - 8))}"
        printf "${CYAN}│${RESET} ${DIM}%s${RESET}%*s${CYAN}│${RESET}\n" \
            "$short_latest" $((W - 4 - ${#short_latest})) ""
    fi

    if [[ -n "$key_line" ]] && [[ "$key_line" != "$latest" ]]; then
        short_key="${key_line:0:$((W - 8))}"
        printf "${CYAN}│${RESET} ${DIM}%s${RESET}%*s${CYAN}│${RESET}\n" \
            "$short_key" $((W - 4 - ${#short_key})) ""
    fi

    printf "${BOLD}${CYAN}├"
    repeat_char "─" $((W - 2))
    printf "┤${RESET}\n"

    printf "${CYAN}│${RESET} ${DIM}Refresh: %ss   Ctrl-C: exit monitor only; IQ-TREE keeps running.${RESET}%*s${CYAN}│${RESET}\n" \
        "$INTERVAL" $((W - 68 - ${#INTERVAL})) ""

    printf "${BOLD}${CYAN}└"
    repeat_char "─" $((W - 2))
    printf "┘${RESET}\n"

    [[ "$status" == "FINISHED" ]] && exit 0
    sleep "$INTERVAL"
done

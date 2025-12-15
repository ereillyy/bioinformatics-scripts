#!/bin/bash
# Summarise per-user CPU and memory usage, sorted by total usage

ps -eo user,%cpu,%mem --no-headers | \
awk '{
    cpu[$1]+=$2; mem[$1]+=$3; total_cpu+=$2; total_mem+=$3
}
END {
    for (u in cpu) {
        sum = cpu[u]/32 + mem[u]
        printf "%-20s %-15.2f %-15.2f %-15.2f\n", u, cpu[u]/32, mem[u], sum
    }
    # Print sum row
    printf "%-20s %-15.2f %-15.2f %-15.2f\n", "TOTAL", total_cpu/32, total_mem, total_cpu/32 + total_mem
}' | sort -k4 -nr | \
awk 'BEGIN {printf "%-20s %-15s %-15s %-15s\n", "USER", "CPU (%)", "MEM (%)", "SUM (CPU+MEM)"} {print}'

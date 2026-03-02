squeue -p epyc -o "%.9u %S %.3C %.16m %.20j %.6k %.2t %.6M %R" | grep Priority | sort -k2 -r | awk '
function mem_to_bytes(v) {
    if (v ~ /G$/) return substr(v,1,length(v)-1) * 1024^3
    if (v ~ /M$/) return substr(v,1,length(v)-1) * 1024^2
    if (v ~ /K$/) return substr(v,1,length(v)-1) * 1024
    return v
}

function bytes_to_human(v) {
    return int(v/1024^3) "G"
}

function reset_block() {
    delete uniq
    unique_count = 0
    sumC = 0
    sumD = 0
    in_block = 0
}

function print_block() {
    if (in_block) {
        printf "%d\t\t%d\t%s\n",
            unique_count,
            sumC,
            bytes_to_human(sumD)
        reset_block()
    }
}

BEGIN {
    reset_block()
}

{
    if ($1 == "reilly01") {
        # end any running collapsed block
        print_block()
        print
    } else {
        # accumulate continuous block
        if (!in_block) in_block = 1

        if (!($1 in uniq)) {
            uniq[$1] = 1
            unique_count++
        }

        sumC += $3
        sumD += mem_to_bytes($4)
    }
}

END {
    print_block()
}
'

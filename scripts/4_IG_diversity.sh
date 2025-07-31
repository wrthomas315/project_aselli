#!/bin/bash

#note run in directory with all files available. 

awk -F',' 'FNR > 1 {print $1","$2}' ML{6..9}_output | sort | uniq > all_unique_CDR3.txt

#file with unique combos
unique_file="all_unique_CDR3.txt"

#name outfiles from ImRep files
files=(ML6_output ML7_output ML8_output ML9_output)

#print the header
echo "CDR3_AA_Seq,Chain_type,ReadCount_ML6,ReadCount_ML7,ReadCount_ML8,ReadCount_ML9" > unique_counts.csv

#process unique file lines except header
tail -n +2 "$unique_file" | while IFS=',' read -r seq chain; do
    #initialize array to hold counts for each file
    counts=()
    for f in "${files[@]}"; do
        #search for exact match
        count=$(awk -F',' -v s="$seq" 'NR>1 && $1 == s {print $3; exit}' "$f")
        if [ -z "$count" ]; then
            count=0
        fi
        counts+=("$count")
    done
    #print new line to output
    echo "$seq,$chain,${counts[0]},${counts[1]},${counts[2]},${counts[3]}" >> unique_counts.csv
done

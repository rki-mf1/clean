#!/bin/bash

FILE=${1}
#LIST=$2

#STAT=$(basename $(basename ${LIST}) .txt)
echo "Filename"$'\t'"Total"$'\t'"mapped_reads"$'\t'"unmapped_reads"$'\t'"proportion_mapped"$'\t'"proportion_unmapped" > ${FILE}.mappedstats.txt


mapped_reads=$(samtools stats ${FILE}.sorted.bam |  grep ^SN | awk 'NR == 7 {print $4}')
unmapped_reads=$(samtools stats ${FILE}.sorted.bam |  grep ^SN | awk 'NR == 9 {print $4}')
total=$(samtools stats ${FILE}.sorted.bam |  grep ^SN | awk 'NR == 1 {print $5}')
proportion_mapped=$(awk "BEGIN {printf \"%.2f\n\", $mapped_reads*100/$total}")
proportion_unmapped=$(awk "BEGIN {printf \"%.2f\n\", $unmapped_reads*100/$total}")
printf "%s\t%s\t%s\t%s\t%s\t%s\n" $FILE $total $mapped_reads $unmapped_reads $proportion_mapped $proportion_unmapped >> ${FILE}.mappedstats.txt
# remove bam file
#rm -rf /MIGE/01_DATA/02_MAPPING/${FILE}*.{bam,bai} 

process bbdukStats {
  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "stats.txt"

  input:
  path bbdukStats

  output:
  path "stats.txt"

  script:
  """
  TOTAL=\$(grep '#Total' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MNUM=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MPER=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$3}')

  FA=\$(awk -F '\\t' '/^[^#]/ {print "\\t\\t"\$2" ("\$3") aligned to "\$1}' ${bbdukStats})

  touch stats.txt
  cat <<EOF >> stats.txt
  \$TOTAL reads; of these:
  \t\$MNUM (\$MPER) were properly paired and mapped; of these:
  \$FA
  EOF
  """
}


process bbduk {
  label 'bbmap'
  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "*.gz" 
  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "log.txt"

  input:
  tuple val(name), file(reads)
  file(db)

  output:
  tuple val(name), file("*clean*.fastq.gz")
  tuple val(name), file("*contamination*.fastq.gz")
  path "stats.txt", emit: stats
  path "log.txt"

  shell:
  """
  # replace the space in the header to retain the full read IDs after mapping (the mapper would split the ID otherwise after the first space)
  # this is working for ENA reads that have at the end of a read id '/1' or '/2'
  EXAMPLE_ID=\$(zcat ${reads[0]} | head -1)
  if [[ \$EXAMPLE_ID == */1 ]]; then 
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      zcat ${reads[0]} | sed 's/ /DECONTAMINATE/g' > ${name}.R1.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${reads[0]} > ${name}.R1.id.fastq
    fi
    if [[ ${reads[1]} =~ \\.gz\$ ]]; then
      zcat ${reads[1]} | sed 's/ /DECONTAMINATE/g' > ${name}.R2.id.fastq
    else
      sed 's/ /DECONTAMINATE/g' ${reads[1]} > ${name}.R2.id.fastq
    fi
  else
    # this is for paried-end SRA reads that don't follow the ENA pattern
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      zcat ${reads[0]} > ${name}.R1.id.fastq
      zcat ${reads[1]} > ${name}.R2.id.fastq
    else
      cp ${reads[0]} ${name}.R1.id.fastq
      cp ${reads[1]} ${name}.R2.id.fastq
    fi
  fi
  
  # bbduk
  echo ${task.memory}
  MEM=\$(echo ${task.memory} | sed 's/ GB//g')
  bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=stats.txt ordered=t k=${params.bbduk_kmer} in=${name}.R1.id.fastq in2=${name}.R2.id.fastq out=${name}.clean.R1.id.fastq out2=${name}.clean.R2.id.fastq outm=${name}.contamination.R1.id.fastq outm2=${name}.contamination.R2.id.fastq

  # restore the original read IDs
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R2.fastq.gz
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R2.fastq.gz
  rm ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq
  
  touch log.txt
  cat <<EOF >> log.txt
Input:\t${reads[0]}, ${reads[1]} 
Host:\t${db}

Clean:\t\t${params.output}/${name}/bbduk/${name}.clean.R1.fastq.gz
\t\t${params.output}/${name}/bbduk/${name}.clean.R2.fastq.gz
Contaminated:\t${params.output}/${name}/bbduk/${name}.contamination.R1.fastq.gz
\t\t${params.output}/${name}/bbduk/${name}.contamination.R2.fastq.gz

# Stay clean!
EOF
"""
}
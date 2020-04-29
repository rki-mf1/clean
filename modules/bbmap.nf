process bbduk {
  label 'bbmap'
  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "*.gz"

  input:
  tuple val(name), file(reads)
  file(db)

  output:
  tuple val(name), file("*clean*.fastq.gz")
  tuple val(name), file("*contamination*.fastq.gz")
  path "bbduk_stats.txt", emit: stats
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
  bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${name}.R1.id.fastq in2=${name}.R2.id.fastq out=${name}.clean.R1.id.fastq out2=${name}.clean.R2.id.fastq outm=${name}.contamination.R1.id.fastq outm2=${name}.contamination.R2.id.fastq

  # restore the original read IDs
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.clean.R2.fastq.gz
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | gzip > ${name}.contamination.R2.fastq.gz
  rm ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq
  """
}
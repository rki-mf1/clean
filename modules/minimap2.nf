/*Comment section: */

process minimap2_fasta {
  label 'minimap2'

  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.gz"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam.bai"

  input: 
    tuple val(name), path(fasta)
    path db

  output:
    val name, emit: name
    path '*.gz'
    path '*.contamination.sorted.bam'
    path '*.contamination.sorted.bam.bai'
    path 'idxstats.tsv', emit: idxstats
    env TOTALREADS, emit: totalreads

  script:
  """
  if [[ ${fasta} =~ \\.gz\$ ]]; then
    TOTALREADS=\$(zgrep '^>' ${fasta} | wc -l)
  else
    TOTALREADS=\$(grep '^>' ${fasta} | wc -l)
  fi

  minimap2 -ax asm5 -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${fasta}

  samtools fasta -f 4 -0 ${name}.clean.fasta ${name}.sam
  samtools fasta -F 4 -0 ${name}.contamination.fasta ${name}.sam
  pigz -p ${task.cpus} ${name}.clean.fasta
  pigz -p ${task.cpus} ${name}.contamination.fasta

  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv

  rm -f ${name}.sam
  """
}

process minimap2_nano {
  label 'minimap2'
  
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.gz"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam.bai"

  input: 
    tuple val(name), path(fastq)
    path db

  output:
    val name, emit: name
    path '*.gz'
    path '*.contamination.sorted.bam'
    path '*.contamination.sorted.bam.bai'
    path 'idxstats.tsv', emit: idxstats
    env TOTALREADS, emit: totalreads

  script:
  """  
  # remove spaces in read IDs to keep them in the later cleaned output
  if [[ ${fastq} =~ \\.gz\$ ]]; then
    zcat ${fastq} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq
    TOTALREADS=\$(zcat ${fastq} | echo \$((`wc -l`/4)))
  else
    sed 's/ /DECONTAMINATE/g' ${fastq} > ${name}.id.fastq
    TOTALREADS=\$(cat ${fastq} | echo \$((`wc -l`/4)))
  fi

  PARAMS="-ax map-ont"
  if [[ ${params.reads_rna} != 'false' ]]; then
    PARAMS="-ax splice -k14"
  fi

  minimap2 \$PARAMS -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq

  samtools fastq -f 4 -0 ${name}.clean.id.fastq ${name}.sam
  samtools fastq -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

  sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | pigz -p ${task.cpus} > ${name}.clean.fastq.gz
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | pigz -p ${task.cpus} > ${name}.contamination.fastq.gz

  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats  ${name}.contamination.sorted.bam > idxstats.tsv

  rm -f ${name}.sam ${name}.clean.id.fastq ${name}.contamination.id.fastq
  """
}

process minimap2_illumina {
  label 'minimap2'

  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.gz"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam.bai"

  input: 
    tuple val(name), path(reads)
    path db
    val mode

  output:
    val name, emit: name
    path '*.gz'
    path '*.contamination.sorted.bam'
    path '*.contamination.sorted.bam.bai'
    path 'idxstats.tsv', emit: idxstats
    env TOTALREADS, emit: totalreads

  script:
  if ( mode == 'paired' ) {
  """
  # replace the space in the header to retain the full read IDs after mapping (the mapper would split the ID otherwise after the first space)
  # this is working for ENA reads that have at the end of a read id '/1' or '/2'
  EXAMPLE_ID=\$(zcat ${reads[0]} | head -1)
  if [[ \$EXAMPLE_ID == */1 ]]; then 
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      zcat ${reads[0]} | sed 's/ /DECONTAMINATE/g' > ${name}.R1.id.fastq
      TOTALREADS_1=\$(zcat ${reads[0]} | echo \$((`wc -l`/4)))
    else
      sed 's/ /DECONTAMINATE/g' ${reads[0]} > ${name}.R1.id.fastq
      TOTALREADS_1=\$(cat ${reads[0]} | echo \$((`wc -l`/4)))
    fi
    if [[ ${reads[1]} =~ \\.gz\$ ]]; then
      zcat ${reads[1]} | sed 's/ /DECONTAMINATE/g' > ${name}.R2.id.fastq
      TOTALREADS_2=\$(zcat ${reads[1]} | echo \$((`wc -l`/4)))
    else
      sed 's/ /DECONTAMINATE/g' ${reads[1]} > ${name}.R2.id.fastq
      TOTALREADS_2=\$(cat ${reads[1]} | echo \$((`wc -l`/4)))
    fi
  else
    # this is for paried-end SRA reads that don't follow the ENA pattern
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      zcat ${reads[0]} > ${name}.R1.id.fastq
      zcat ${reads[1]} > ${name}.R2.id.fastq
      TOTALREADS_1=\$(zcat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALREADS_2=\$(zcat ${reads[1]} | echo \$((`wc -l`/4)))
    else
      mv ${reads[0]} ${name}.R1.id.fastq
      mv ${reads[1]} ${name}.R2.id.fastq
      TOTALREADS_1=\$(cat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALREADS_2=\$(cat ${reads[1]} | echo \$((`wc -l`/4)))
    fi
  fi
  TOTALREADS=\$(( TOTALREADS_1+TOTALREADS_2 ))

  # Use samtools -F 2 to discard only reads mapped in proper pair:
  minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${name}.R1.id.fastq ${name}.R2.id.fastq
  
  samtools fastq -F 2 -1 ${name}.clean.R1.id.fastq -2 ${name}.clean.R2.id.fastq ${name}.sam
  samtools fastq -f 2 -1 ${name}.contamination.R1.id.fastq -2 ${name}.contamination.R2.id.fastq ${name}.sam

  samtools view -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv

  # restore the original read IDs
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | pigz -p ${task.cpus} > ${name}.clean.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.clean.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | pigz -p ${task.cpus} > ${name}.clean.R2.fastq.gz
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R1.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/1"}else{print \$0};LINE++;}' | pigz -p ${task.cpus} > ${name}.contamination.R1.fastq.gz 
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.R2.id.fastq | awk 'BEGIN{LINE=0};{if(LINE % 4 == 0 || LINE == 0){print \$0"/2"}else{print \$0};LINE++;}' | pigz -p ${task.cpus} > ${name}.contamination.R2.fastq.gz

  # remove intermediate files
  rm -f ${name}.R1.id.fastq ${name}.R2.id.fastq ${name}.clean.R1.id.fastq ${name}.clean.R2.id.fastq ${name}.contamination.R1.id.fastq ${name}.contamination.R2.id.fastq ${name}.sam
  """
  } else {
  """
  # remove spaces in read IDs to keep them in the later cleaned output
  if [[ ${reads} =~ \\.gz\$ ]]; then
    zcat ${reads} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq
    TOTALREADS=\$(zcat ${reads} | echo \$((`wc -l`/4)))
  else
    sed 's/ /DECONTAMINATE/g' ${reads} > ${name}.id.fastq
    TOTALREADS=\$(cat ${reads} | echo \$((`wc -l`/4)))
  fi

  minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq
  
  samtools fastq -f 4 -0 ${name}.clean.id.fastq ${name}.sam
  samtools fastq -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv

  sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | pigz -p ${task.cpus} > ${name}.clean.fastq.gz
  sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | pigz -p ${task.cpus} > ${name}.contamination.fastq.gz

  # remove intermediate files
  rm -f ${name}.id.fastq ${name}.clean.id.fastq ${name}.contamination.id.fastq ${name}.sam
  """
  }
}

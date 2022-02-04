process bbduk {
  label 'bbmap'
  
  input:
  tuple val(name), path(reads)
  path db

  output:
  val name, emit: name
  tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
  tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads
  tuple val(name), path('bbduk_stats.txt'), emit: stats

  script:
  if ( params.mode == 'paired' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads[0]} in2=${reads[1]} out=${reads[0].baseName}.clean.fastq out2=${reads[1].baseName}.clean.fastq outm=${reads[0].baseName}.contamination.fastq outm2=${reads[1].baseName}.contamination.fastq
    """
  } else if ( params.mode == 'single' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads} out=${reads.baseName}.clean.fastq outm=${reads.baseName}.contamination.fastq
    """
  } else {
    error "Invalid mode: ${params.mode}"
  }
  stub:
  """
  touch bbduk_stats.txt
  touch ${reads.baseName}.contamination.fastq ${reads.baseName}.clean.fastq
  """
}

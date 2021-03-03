process bbduk {
  label 'bbmap'
  
  input:
  tuple val(name), path(reads)
  path db
  val mode

  output:
  val name, emit: name
  tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
  tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads
  path 'bbduk_stats.txt', emit: stats

  script:
  if ( mode == 'paired' ) {
  """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads[0]} in2=${reads[1]} out=${reads[0].baseName}.clean.fastq out2=${reads[1].baseName}.clean.fastq outm=${reads[0].baseName}.contamination.fastq outm2=${reads[1].baseName}.contamination.fastq
  """
  } else {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads} out=${reads.baseName}.clean.fastq outm=${reads.baseName}.contamination.fastq
    """
  }
}

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
  if ( params.lib_pairedness == 'paired' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads[0]} in2=${reads[1]} out=${reads[0].baseName}.clean.fastq out2=${reads[1].baseName}.clean.fastq outm=${reads[0].baseName}.contamination.fastq outm2=${reads[1].baseName}.contamination.fastq
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads} out=${name}.clean.fastq outm=${name}.contamination.fastq
    """
  } else {
    error "Invalid mode: ${params.lib_pairedness}"
  }
  stub:
  if ( params.lib_pairedness == 'paired' ) {
    """
    touch bbduk_stats.txt
    touch ${reads[0].baseName}.clean.fastq ${reads[0].baseName}.contamination.fastq ${reads[1].baseName}.clean.fastq ${reads[1].baseName}.contamination.fastq
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    touch bbduk_stats.txt
    touch ${name}.contamination.fastq ${name}.clean.fastq
    """
  }
}

process bbduk {
  label 'bbmap'
  
  publishDir "${params.output}/${params.tool}/${name}", mode: 'copy', pattern: "*.gz", enabled: false
  
  input:
  tuple val(name), path(reads)
  path db

  output:
  val name, emit: name
  tuple val(name), val('clean'), path('*clean.fastq.gz'), emit: cleaned_reads
  tuple val(name), val('contamination'), path('*contamination.fastq.gz'), emit: contaminated_reads
  tuple val(name), path("${name}_bbduk_stats.txt"), emit: stats

  script:
  if ( params.lib_pairedness == 'paired' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=${name}_bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads[0]} in2=${reads[1]} out=${reads[0].baseName}.clean.fastq out2=${reads[1].baseName}.clean.fastq outm=${reads[0].baseName}.contamination.fastq outm2=${reads[1].baseName}.contamination.fastq

    gzip --no-name *.clean.fastq
    gzip --no-name *.contamination.fastq
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=${name}_bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads} out=${name}.clean.fastq outm=${name}.contamination.fastq

    gzip --no-name *.clean.fastq
    gzip --no-name *.contamination.fastq
    """
  } else {
    error "Invalid mode: ${params.lib_pairedness}"
  }
  stub:
  if ( params.lib_pairedness == 'paired' ) {
    """
    touch ${name}_bbduk_stats.txt
    touch ${reads[0].baseName}.clean.fastq.gz ${reads[0].baseName}.contamination.fastq.gz ${reads[1].baseName}.clean.fastq.gz ${reads[1].baseName}.contamination.fastq.gz
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    touch ${name}_bbduk_stats.txt
    touch ${name}.contamination.fastq.gz ${name}.clean.fastq.gz
    """
  }
}

process bbduk {
  label 'bbmap'

  publishDir (
    path: "${params.output}/intermediate/${map_target}",
    mode: params.publish_dir_mode,
    pattern: "*.{clean,contamination}.fastq.gz",
    enabled: !params.no_intermediate,
    overwrite: false
  )

  // When using `--keep`, we need to do further processing before we have
  // the final clean and removed data sets.
  if ( !params.keep ) {
    publishDir (
      path: params.output,
      overwrite: false,
      mode: params.publish_dir_mode,
      pattern: "*.gz",
      saveAs: { fn ->
            fn.endsWith('.clean.fastq.gz') ? "clean/${fn}".replaceAll(~'.fastq.clean.fastq.gz$', '.fastq.gz') :
            fn.endsWith('.contamination.fastq.gz') ? "removed/${fn}".replaceAll(~'.fastq.contamination.fastq.gz$', '.fastq.gz') :
            fn
      }
    )
  }

  input:
  tuple val(name), path(reads)
  path db
  val(map_target)

  output:
  val name, emit: name
  tuple val(name), val('unmapped'), path('*clean.fastq.gz'), emit: cleaned_reads
  tuple val(name), val('mapped'), path('*contamination.fastq.gz'), emit: contaminated_reads
  tuple val(name), path("${name}.bbduk_stats.txt"), emit: stats

  script:
  if ( params.lib_pairedness == 'paired' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=${name}.bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads[0]} in2=${reads[1]} out=${reads[0].baseName}.clean.fastq out2=${reads[1].baseName}.clean.fastq outm=${reads[0].baseName}.contamination.fastq outm2=${reads[1].baseName}.contamination.fastq

    gzip --no-name *.clean.fastq
    gzip --no-name *.contamination.fastq
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    MEM=\$(echo ${task.memory} | sed 's/ GB//g')
    bbduk.sh -Xmx\${MEM}g ref=${db} threads=${task.cpus} stats=${name}.bbduk_stats.txt ordered=t k=${params.bbduk_kmer} in=${reads} out=${name}.clean.fastq outm=${name}.contamination.fastq

    gzip --no-name *.clean.fastq
    gzip --no-name *.contamination.fastq
    """
  } else {
    error "Invalid mode: ${params.lib_pairedness}"
  }
  stub:
  if ( params.lib_pairedness == 'paired' ) {
    """
    touch ${name}.bbduk_stats.txt
    touch ${reads[0].baseName}.clean.fastq.gz ${reads[0].baseName}.contamination.fastq.gz ${reads[1].baseName}.clean.fastq.gz ${reads[1].baseName}.contamination.fastq.gz
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    touch ${name}.bbduk_stats.txt
    touch ${name}.contamination.fastq.gz ${name}.clean.fastq.gz
    """
  }
}

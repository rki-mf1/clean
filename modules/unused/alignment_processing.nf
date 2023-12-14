process filter_un_mapped_alignments {
  label 'minimap2'

  input:
    tuple val(name), path(sam), path(reads)

  output:
    tuple val(name), val('clean'), path('*clean.fast{q,a}'), emit: cleaned_reads
    tuple val(name), val('mapped'), path('*mapped.fast{q,a}'), emit: contaminated_reads

  script:
  if ( params.mode == 'paired' ) {
    """
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -@ ${task.cpus} -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -@ ${task.cpus} -f 2 -1 ${reads[0].baseName}.mapped.fastq -2 ${reads[1].baseName}.mapped.fastq ${name}.sam
    """
  } else if ( params.mode == 'single' ) {
    dtype = (params.input_type == 'fasta') ? 'a' : 'q'
    """
    samtools fast${dtype} -@ ${task.cpus} -f 4 -0 ${reads.baseName}.clean.fast${dtype} ${sam}
    samtools fast${dtype} -@ ${task.cpus} -F 4 -0 ${reads.baseName}.mapped.fast${dtype} ${sam}
    """
  } else {
    error "Invalid mode: ${params.mode}"
  }
  stub:
  """
  touch ${reads.baseName}.clean.fasta ${reads.baseName}.mapped.fasta ${reads.baseName}.clean.fastq ${reads.baseName}.mapped.fastq
  """
}


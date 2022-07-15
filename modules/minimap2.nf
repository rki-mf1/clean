process minimap2 {
  label 'minimap2'

  input: 
    tuple val(name), path(input)
    path (db)

  output:
    tuple val(name), val('raw'), path("${name}.bam"), emit: bam // input just for naming

  script:
  // -N is an internal algorithm option. It controls how many candidates alignment to extend. --secondary is an output option.
  if ( params.input_type == 'nano' ) {
    params = params.reads_rna ? "-ax splice -k14" : "-ax map-ont"
    """
    minimap2 ${params} -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} ${db} ${input} | samtools view -bhS -@ ${task.cpus} > ${name}.bam
    """
  } else if ( params.input_type.contains('illumina') ) {
    """
    minimap2 -ax sr -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} ${db} ${input} | samtools view -bhS -@ ${task.cpus} > ${name}.bam
    """
  } else if ( params.input_type == 'fasta' ){
    """
    minimap2 -ax asm5 -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} ${db} ${input} | samtools view -bhS -@ ${task.cpus} > ${name}.bam
    """
  } else {
    error "Unknown input_type: ${params.input_type}"
  }
  stub:
  """
  touch ${name}.bam
  """
}

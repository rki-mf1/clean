// -k and -w are baked into the index, so the index has to be built with the
// same preset that is later used for the mapping
def minimap2_preset() {
  if ( params.input_type == 'nano' ) {
    return params.reads_rna ? '-x splice -k14' : '-x map-ont'
  } else if ( params.input_type == 'pacbio' ) {
    return params.reads_rna ? '-x splice -k14' : '-x map-pb'
  } else if ( params.input_type.contains('illumina') ) {
    return '-x sr'
  } else if ( params.input_type == 'fasta' ) {
    return '-x asm5'
  } else {
    error "Unknown input_type: ${params.input_type}"
  }
}

process minimap2_index {
  label 'minimap2'

  input:
    path (fasta)

  output:
    path 'db.mmi'

  script:
  // without this every sample rebuilds the index, which dominates the runtime
  // for large hosts
  """
  minimap2 ${minimap2_preset()} -t ${task.cpus} -d db.mmi ${fasta}
  """
  stub:
  """
  touch db.mmi
  """
}

process minimap2 {
  label 'minimap2'

  input:
    tuple val(name), path(input)
    path (db)

  output:
    tuple val(name), val('raw'), path("${name}.bam"), emit: bam // input just for naming

  script:
  // -N is an internal algorithm option. It controls how many candidates alignment to extend. --secondary is an output option.
  // --split-prefix keeps the output correct when the index is split into parts
  """
  minimap2 -a ${minimap2_preset()} -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} ${db} ${input} | samtools view -bhS -@ ${task.cpus} > ${name}.bam
  """
  stub:
  """
  touch ${name}.bam
  """
}

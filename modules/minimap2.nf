process minimap2 {
  label 'minimap2'

  input: 
    tuple val(name), path(input)
    path (db)

  output:
    tuple val(name), env(TOTALCONTIGS), emit: num_contigs optional true
    tuple val(name), path('*.sam'), path(input), emit: sam // input just for naming

  script:
  if ( params.seq_type == 'nano' ) {
    """
    PARAMS="-ax map-ont"
    if [[ ${params.reads_rna} != 'false' ]]; then
      PARAMS="-ax splice -k14"
    fi
    
    minimap2 \$PARAMS -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${input}
    """
  } else if ( params.seq_type == 'illumina' ) {
    if ( params.mode == 'paired' ) {
      """
      minimap2 -ax sr -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${input[0]} ${input[1]}
      """
    } else {
      """
      minimap2 -ax sr -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${input}
      """
    }
  } else if ( params.seq_type == 'fasta' ){
      """
      if [[ ${fasta} =~ \\.gz\$ ]]; then
        TOTALCONTIGS=\$(zgrep '^>' ${fasta} | wc -l)
      else
        TOTALCONTIGS=\$(grep '^>' ${fasta} | wc -l)
      fi

      minimap2 -ax asm5 -N 5 --split-prefix tmp --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${fasta}
      """
  } else {
    error "Unknown seq_type: ${params.seq_type}"
  }
  stub:
  """
  TOTALCONTIGS=42
  touch ${name}.sam
  """
}

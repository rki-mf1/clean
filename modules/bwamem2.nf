process bwamem2_index {
  label 'bwamem2_index'

  input:
    path(fasta)

  output:
    path('bwamem2') , emit: index

  script:
  """
  mkdir bwamem2
  bwa-mem2 \\
    index \\
    -p bwamem2/db \\
    $fasta
  """

  stub:
  """
  mkdir bwamem2

  touch bwamem2/db.{0123,amb,ann,bwt.2bit.64,pac}
  """
}

process bwamem2 {
  label 'bwamem2'

  input:
  tuple val(name), path(input)
  path(db_index)
  path(db)


  output:
  tuple val(name), val('raw'), path("${name}.bam"), emit: bam // input just for naming

  script:
  """
  INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
  bwa-mem2 mem \\
    -t $task.cpus \\
    \$INDEX \\
    $input \\
    | samtools view -bhS -@ $task.cpus > ${name}.bam
  """
  stub:
  """
  touch ${name}.bam
  """
}

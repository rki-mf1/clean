process bwa_index {
  label 'bwa'
  
  input:
    path(fasta)
  
  output:
    path(bwa) , emit: index
  
  script:
  """
  mkdir bwa
  bwa \\
    index \\
    -p bwa/db \\
    $fasta
  """
  
  stub:
  """
  mkdir bwa
  
  touch bwa/db.{amb,ann,bwt,pac,sa}
  """
}

process bwa {
  label 'bwa'

  input: 
  tuple val(name), path(input)
  path(db_index)
  path(db)


  output:
  tuple val(name), val('raw'), path("${name}.bam"), emit: bam // input just for naming

  script:
  """
  INDEX=`find -L ./ -name "*.amb" | sed 's/\\.amb\$//'`
  bwa mem \\
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

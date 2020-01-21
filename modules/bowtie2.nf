/*Comment section: */

process bowtie2_illumina {
  label 'bowtie2'
  publishDir "${params.output}/${name}/bowtie_ill_martin_extraction", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(reads)
    file(genome)
    file(index)

  output:
    file("*.gz")

  script:
    """
    NAME=${name}
    if [ "${params.phix}" != "false" ]; then
      NAME=${name}_phix
    fi

    bowtie2 -x \${NAME} -1 ${reads[0]} -2 ${reads[1]} -S ${name}.sam

    samtools fastq -F 2 -1 ${name}.clean.R1.fastq -2 ${name}.clean.R2.fastq ${name}.sam

    rm *.gz
    gzip ${name}.clean.R1.fastq
    gzip ${name}.clean.R2.fastq

    """
}

process bowtie2_illumina_ebi_extraction {
  label 'bowtie2'
  publishDir "${params.output}/${name}/bowtie_ill_ebi_extraction", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(reads)
    file(db)

  output:
    file("*.gz")

  script:
    """
    NAME=${name}
    if [ "${params.phix}" != "false" ]; then
      NAME=${name}_phix
    fi

    bowtie2 -x \${NAME} -1 ${reads[0]} -2 ${reads[1]} -S ${name}.sam

    samtools fastq -f 12 -F 256 -1 ${name}.clean.R1.fastq -2 ${name}.clean.R2.fastq ${name}.sam

    rm *.gz
    gzip ${name}.clean.R1.fastq
    gzip ${name}.clean.R2.fastq

    """
}

/*Comment section: */

process minimap2_fasta {
  label 'minimap2'
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(fasta)
    file(db)

  output:
    file("*.gz")

  script:
    """
    minimap2 -ax asm5 -t ${task.cpus} -o ${name}.sam ${db} ${fasta}
    samtools fasta -f 4 -0 ${name}.clean.fasta ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.fasta ${name}.sam
    gzip -f ${name}.clean.fasta
    gzip -f ${name}.contamination.fasta
    rm ${name}.sam
    """
}

process minimap2_nano {
  label 'minimap2'
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "*.gz" 

  input: 
    tuple val(name), file(fastq)
    file(db)

  output:
    file("*.gz")

  script:
    """
    # remove spaces in read IDs to keep them in the later cleaned output
    zcat ${fastq} | sed 's/ /DECONTAMINATE/g' > ${name}.id.fastq

    minimap2 -ax map-ont -t ${task.cpus} -o ${name}.sam ${db} ${name}.id.fastq
    samtools fasta -f 4 -0 ${name}.clean.id.fastq ${name}.sam
    samtools fasta -F 4 -0 ${name}.contamination.id.fastq ${name}.sam

    sed 's/DECONTAMINATE/ /g' ${name}.clean.id.fastq | gzip > ${name}.clean.fastq.gz
    sed 's/DECONTAMINATE/ /g' ${name}.contamination.id.fastq | gzip > ${name}.contamination.fastq.gz
     
    rm ${name}.sam ${name}.clean.id.fastq ${name}.contamination.id.fastq ${name}.id.fastq
    """
}

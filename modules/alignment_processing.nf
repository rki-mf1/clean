process filter_un_mapped_alignments {
  label 'minimap2'

  input:
    tuple val(name), path(sam), path(reads)
    val(mode)

  output:
    tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
    tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads

  script:
  if ( mode == 'paired' ) {
    """
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -f 2 -1 ${reads[0].baseName}.contamination.fastq -2 ${reads[1].baseName}.contamination.fastq ${name}.sam
    """
  } else {
    """
    samtools fastq -f 4 -0 ${reads.baseName}.clean.fastq ${sam}
    samtools fastq -F 4 -0 ${reads.baseName}.contamination.fastq ${sam}
    """
  }
}

process make_contamination_bam {
  label 'minimap2'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "*.contamination.sorted.bam*"

  input:
    tuple val(name), path(sam), path(reads)
    val(mode)
    val(tool)

  output:
    tuple val(name), path ('*.contamination.sorted.bam'), emit: contamination_bam
    tuple val(name), path ('*.contamination.sorted.bam.bai'), emit: contamination_bai
    tuple val(name), path ('idxstats.tsv'), emit: idxstats

  script:
  if ( mode == 'paired' ) {
    """
    samtools view -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv
    """
  } else {
    """
    samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats  ${name}.contamination.sorted.bam > idxstats.tsv
    """
  }
}

process filter_soft_clipped_alignments {
  label 'teloclip'

  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  path (db)
  path (fai)
  val (max)
  val (tool)
  
  output:
  path ('*.ambiguous.bam')
  path ('*.ambiguous.bam.bai')
  
  script:
  """
  unpigz -f -k -p ${task.cpus} ${db}
  samtools view -h ${bam} | teloclip --minClip ${max} --ref ${fai} | samtools sort > ${name}.ambiguous.bam
  samtools index ${name}.ambiguous.bam
  rm -f ${db.baseName}
  """
}
process filter_un_mapped_alignments {
  label 'minimap2'

  input:
    tuple val(name), path(sam), path(reads)
    val(mode)

  output:
    tuple val(name), val('clean'), path('*clean.fast{q,a}'), emit: cleaned_reads
    tuple val(name), val('mapped'), path('*mapped.fast{q,a}'), emit: contaminated_reads

  script:
  if ( mode == 'paired' ) {
    """
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -f 2 -1 ${reads[0].baseName}.mapped.fastq -2 ${reads[1].baseName}.mapped.fastq ${name}.sam
    """
  } else if ( mode == 'single' || mode == 'fasta' ) {
    dtype = (mode == 'single') ? 'q' : 'a'
    """
    samtools fast${dtype} -f 4 -0 ${reads.baseName}.clean.fast${dtype} ${sam}
    samtools fast${dtype} -F 4 -0 ${reads.baseName}.mapped.fast${dtype} ${sam}
    """
  } else {
    error "Invalid mode: ${mode}"
  }
}

process make_mapped_bam {
  label 'minimap2'

  publishDir "${params.output}/${name}/${params.tool}", mode: 'copy', pattern: "*.mapped.bam*"

  input:
    tuple val(name), path(sam), path(reads)

  output:
    tuple val(name), path ('*.mapped.bam'), emit: contamination_bam
    tuple val(name), path ('*.mapped.bam.bai'), emit: contamination_bai
    tuple val(name), path ('idxstats.tsv'), emit: idxstats

  script:
  if ( params.mode == 'paired' ) {
    """
    samtools view -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.mapped.bam --threads ${task.cpus}
    samtools index ${name}.mapped.bam
    samtools idxstats ${name}.mapped.bam > idxstats.tsv
    """
  } else {
    """
    samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.mapped.bam --threads ${task.cpus}
    samtools index ${name}.mapped.bam
    samtools idxstats  ${name}.mapped.bam > idxstats.tsv
    """
  }
}

process filter_soft_clipped_alignments {
  label 'samclipy'
  label 'smallTask'

  publishDir "${params.output}/${name}/${params.tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  val (minClip)
  
  output:
  tuple val(name), val('ambiguous'), path ('*.ambiguous.bam'), emit: bam_am
  tuple val(name), val('contamination'), path ('*.contamination.bam'), emit: bam_unam
  tuple val(name), path ('*.bam.bai')
  
  script:
  """
  git clone https://github.com/MarieLataretu/samclipy.git --branch v0.0.2 || git clone git@github.com:MarieLataretu/samclipy.git --branch v0.0.2 
  samtools view -h ${bam} | python samclipy/samclipy.py --invert --minClip ${minClip} | samtools sort > ${name}.ambiguous.bam
  samtools view -h ${bam} | python samclipy/samclipy.py --minClip ${minClip} | samtools sort > ${name}.contamination.bam
  samtools index ${name}.ambiguous.bam
  samtools index ${name}.contamination.bam
  """
}

process fastq_from_bam {
  label 'minimap2'

  input:
  tuple val(name), val(type), path(bam)

  output:
  tuple val(name), val(type), path('*.fastq')

  script:
  if ( params.mode == 'paired' ) {
    """
    samtools fastq -@ ${task.cpus} -1 ${bam.baseName}_1.fastq -2 ${bam.baseName}_2.fastq -s ${bam.baseName}_singleton.fastq ${bam}
    """
  } else {
    """
    samtools fastq -@ ${task.cpus} -0 ${bam.baseName}.fastq ${bam}
    """
  }
}

process idxstats_from_bam {
  label 'minimap2'

  input:
  tuple val(name), path(bam)

  output:
  path('*_idxstats.tsv')

  script:
  """
  samtools idxstats ${bam} > ${bam.baseName}_idxstats.tsv
  """
}
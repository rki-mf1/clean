process split_bam {
  label 'minimap2'

  input:
    tuple val(name), path(bam)

  output:
    tuple val(name), val('mapped'), path("${name}.mapped.bam"), emit: mapped
    tuple val(name), val('unmapped'), path("${name}.unmapped.bam"), emit: unmapped

  script:
  // includes supplementary alignments included (chimeric alignment, sometimes also not linear )
  if ( params.lib_pairedness == 'paired' ){
    """
    samtools view -h -@ ${task.cpus} -f 2 ${bam} | samtools sort -o ${name}.mapped.bam -@ ${task.cpus}
    samtools view -h -@ ${task.cpus} -F 2 ${bam} | samtools sort -o ${name}.unmapped.bam -@ ${task.cpus}
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    samtools view -h -@ ${task.cpus} -F 4 ${bam} | samtools sort -o ${name}.mapped.bam -@ ${task.cpus}
    samtools view -h -@ ${task.cpus} -f 4 ${bam} | samtools sort -o ${name}.unmapped.bam -@ ${task.cpus}
    """
  } else { error "Invalid pairedness: ${params.lib_pairedness}" }
  stub:
  """
  touch ${name}.mapped.bam ${name}.unmapped.bam
  """
}

process merge_bam {
  label 'minimap2'

  input:
    tuple val(name), val(type), path(bam)

  output:
   tuple val(name), val(type), path("${bam[0].baseName}_merged.bam")
  
  script:
  """
  samtools merge -@ ${task.cpus} ${bam[0].baseName}_merged.bam ${bam} # first bam is output
  """
  stub:
  """
  touch ${bam[0].baseName}_merged.bam
  """
}

process filter_soft_clipped_alignments {
  label 'samclipy'
  label 'smallTask'

  publishDir "${params.output}/${params.tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  val (minClip)
  
  output:
  tuple val(name), val('unmapped'), path ('*.softclipped.bam'), emit: bam_clipped
  tuple val(name), val('mapped'), path ('*.passedclipped.bam'), emit: bam_ok_clipped
  tuple val(name), path ('*.bam.bai')
  
  script:
  """
  git clone https://github.com/MarieLataretu/samclipy.git --branch v0.0.2 || git clone git@github.com:MarieLataretu/samclipy.git --branch v0.0.2 
  samtools view -h ${bam} | python samclipy/samclipy.py --invert --minClip ${minClip} | samtools sort > ${name}.softclipped.bam
  samtools view -h ${bam} | python samclipy/samclipy.py --minClip ${minClip} | samtools sort > ${name}.passedclipped.bam
  samtools index ${name}.softclipped.bam
  samtools index ${name}.passedclipped.bam
  """
  stub:
  """
  touch ${name}.softclipped.bam ${name}.passedclipped.bam ${name}.softclipped.bam.bai ${name}.passedclipped.bam.bai
  """
}

process filter_true_dcs_alignments {
  label 'bed_samtools'

  publishDir "${params.output}/${params.tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  path (dcs_ends_bed)

  output:
  tuple val(name), val('mapped'), path ("${name}_no_dcs.bam"), emit: no_dcs
  tuple val(name), val('mapped'), path ("${name}_true_dcs.bam"), emit: true_dcs
  tuple val(name), val('unmapped'), path ("${name}_false_dcs.bam"), emit: false_dcs
  tuple val(name), path ('*.bam.bai')
  path('dcs.bam')

  script:
  """
  # true spike in: 1-65 || 1-92; 3513-3560 (len 48)
  samtools view -b -h -e 'rname=="Lambda_3.6kb"' ${bam} > dcs.bam
  samtools view -b -h -e 'rname!="Lambda_3.6kb"' ${bam} > ${name}_no_dcs.bam
  bedtools intersect -wa -ubam -a dcs.bam -b ${dcs_ends_bed} > ${name}_true_dcs.bam
  bedtools intersect -v -ubam -a dcs.bam -b ${dcs_ends_bed} > ${name}_false_dcs.bam
  samtools index dcs.bam
  samtools index ${name}_no_dcs.bam
  samtools index ${name}_true_dcs.bam
  samtools index ${name}_false_dcs.bam
  """ 
  stub:
  """
  touch ${name}_no_dcs.bam ${name}_true_dcs.bam ${name}_false_dcs.bam ${name}_no_dcs.bam.bai ${name}_true_dcs.bam.bai ${name}_false_dcs.bam.bai
  """
}

process fastq_from_bam {
  label 'minimap2'

  input:
  tuple val(name), val(type), path(bam)

  output:
  tuple val(name), val(type), path('*.fastq')

  script:
  if ( params.lib_pairedness == 'paired' ) {
    """
    samtools fastq -@ ${task.cpus} -1 ${bam.baseName}_1.fastq -2 ${bam.baseName}_2.fastq -s ${bam.baseName}_singleton.fastq ${bam}
    """
  } else if ( params.lib_pairedness == 'single' ) {
    """
    samtools fastq -@ ${task.cpus} -0 ${bam.baseName}.fastq ${bam}
    """
  } else {
    error "Invalid pairedness: ${params.lib_pairedness}"
  }
  stub:
  """
  touch ${bam.baseName}_1.fastq ${bam.baseName}_2.fastq
  """
}

process idxstats_from_bam {
  label 'minimap2'

  input:
  tuple val(name), path(bam), path(bai)

  output:
  tuple val(name), path('*_idxstats.tsv')

  script:
  """
  samtools idxstats ${bam} > ${bam.baseName}_idxstats.tsv
  """
  stub:
  """
  touch ${bam.baseName}_idxstats.tsv
  """
}

process flagstats_from_bam {
  label 'minimap2'

  publishDir "${params.output}/minimap2", mode: 'copy', pattern: "${bam.baseName}_flagstats.txt" 

  input:
  tuple val(name), path(bam), path(bai)

  output:
  tuple val(name), path('*_flagstats.txt')

  script:
  """
  samtools flagstats ${bam} > ${bam.baseName}_flagstats.txt
  """
  stub:
  """
  touch ${bam.baseName}_flagstats.txt
  """
}

process sort_bam {
  label 'minimap2'

  input:
  tuple val(name), path(bam)

  output:
  tuple val(name), path("${bam.baseName}_sorted.bam")

  script:
  """
  samtools sort -@ ${task.cpus} ${bam} > ${bam.baseName}_sorted.bam
  """
  stub:
  """
  touch ${bam.baseName}_sorted.bam
  """
}

process index_bam {
  label 'minimap2'

  input:
  tuple val(name), path(bam)

  output:
  tuple val(name), path(bam), path('*.bai')

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  """
  stub:
  """
  touch ${bam}.bai
  """
}
process filter_un_mapped_alignments {
  label 'minimap2'

  input:
    tuple val(name), path(sam), path(reads)

  output:
    tuple val(name), val('clean'), path('*clean.fast{q,a}'), emit: cleaned_reads
    tuple val(name), val('mapped'), path('*mapped.fast{q,a}'), emit: contaminated_reads

  script:
  if ( params.mode == 'paired' ) {
    """
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -@ ${task.cpus} -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -@ ${task.cpus} -f 2 -1 ${reads[0].baseName}.mapped.fastq -2 ${reads[1].baseName}.mapped.fastq ${name}.sam
    """
  } else if ( params.mode == 'single' ) {
    dtype = (params.seq_type == 'fasta') ? 'a' : 'q'
    """
    samtools fast${dtype} -@ ${task.cpus} -f 4 -0 ${reads.baseName}.clean.fast${dtype} ${sam}
    samtools fast${dtype} -@ ${task.cpus} -F 4 -0 ${reads.baseName}.mapped.fast${dtype} ${sam}
    """
  } else {
    error "Invalid mode: ${params.mode}"
  }
  stub:
  """
  touch ${reads.baseName}.clean.fasta ${reads.baseName}.mapped.fasta ${reads.baseName}.clean.fastq ${reads.baseName}.mapped.fastq
  """
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
    samtools view -@ ${task.cpus} -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.mapped.bam --threads ${task.cpus}
    samtools index ${name}.mapped.bam
    samtools idxstats ${name}.mapped.bam > idxstats.tsv
    """
  } else if ( params.mode == 'single' ) {
    """
    samtools view -@ ${task.cpus} -b -F 2052 ${name}.sam | samtools sort -o ${name}.mapped.bam --threads ${task.cpus}
    samtools index ${name}.mapped.bam
    samtools idxstats  ${name}.mapped.bam > idxstats.tsv
    """
  } else {
    error "Invalid mode: ${params.mode}"
  }
  stub:
  """
  touch ${name}.mapped.bam ${name}.mapped.bam.bai idxstats.tsv
  """
}

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
    path(bam)

  output:
    path("merged.bam")
  
  script:
  """
  samtools merge -@ ${task.cpus} -o merged.bam ${bam} # -h FILE ?
  """
  stub:
  """
  touch merged.bam
  """
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
  stub:
  """
  touch ${name}.ambiguous.bam ${name}.contamination.bam ${name}.ambiguous.bam.bai ${name}.contamination.bam.bai
  """
}

process filter_true_dcs_alignments {
  label 'bed_samtools'

  publishDir "${params.output}/${name}/${params.tool}", mode: 'copy', pattern: "*.bam*"

  input:
  tuple val(name), path (bam)
  path (dcs_ends_bed)

  output:
  tuple val(name), path ("${name}_no_dcs.bam"), emit: no_dcs
  tuple val(name), path ("${name}_true_dcs.bam"), emit: true_dcs
  tuple val(name), path ("${name}_false_dcs.bam"), emit: false_dcs

  script:
  """
  # true spike in: 1-65 || 1-92; 3513-3560 (len 48)
  samtools view -b -h -e 'rname=="Lambda_3.6kb"' ${bam} > tmp.bam
  samtools view -b -h -e 'rname!="Lambda_3.6kb"' ${bam} > _no_dcs.bam
  bedtools intersect -wa -ubam -header -a tmp.bam -b ${dcs_ends_bed} > ${name}_true_dcs.bam
  bedtools intersect -v -ubam -header -a tmp.bam -b ${dcs_ends_bed} > ${name}_false_dcs.bam
  """ 
  stub:
  """
  touch ${name}_no_dcs.bam ${name}_true_dcs.bam ${name}_false_dcs.bam
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
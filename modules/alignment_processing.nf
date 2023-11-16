process split_bam {
  label 'minimap2'

  input:
    tuple val(name), val(type), path(bam)

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

  publishDir "${params.output}/intermediate/soft-clipped", mode: params.publish_dir_mode, pattern: "*.bam*", overwrite: false

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

  publishDir "${params.output}/intermediate/true-dcs", mode: params.publish_dir_mode, pattern: "*.bam*", overwrite: false

  input:
  tuple val(name), path (bam)
  path (dcs_ends_bed)

  output:
  tuple val(name), val('mapped'), path ("${name}.no_dcs.bam"), emit: no_dcs
  tuple val(name), val('mapped'), path ("${name}.true_dcs.bam"), emit: true_dcs
  tuple val(name), val('unmapped'), path ("${name}.false_dcs.bam"), emit: false_dcs
  tuple val(name), path ('*.bam.bai')
  path('dcs.bam')

  script:
  """
  # true spike in: 1-65 || 1-92; 3513-3560 (len 48)
  samtools view -b -h -e 'rname=="Lambda_3.6kb"' ${bam} > dcs.bam
  samtools view -b -h -e 'rname!="Lambda_3.6kb"' ${bam} > ${name}.no_dcs.bam
  bedtools intersect -wa -ubam -a dcs.bam -b ${dcs_ends_bed} > ${name}.true_dcs.bam
  bedtools intersect -v -ubam -a dcs.bam -b ${dcs_ends_bed} > ${name}.false_dcs.bam
  samtools index dcs.bam
  samtools index ${name}.no_dcs.bam
  samtools index ${name}.true_dcs.bam
  samtools index ${name}.false_dcs.bam
  """ 
  stub:
  """
  touch ${name}.no_dcs.bam ${name}.true_dcs.bam ${name}.false_dcs.bam ${name}.no_dcs.bam.bai ${name}.true_dcs.bam.bai ${name}.false_dcs.bam.bai
  """
}

process fastq_from_bam {
  label 'minimap2'

  publishDir (
    path: "${params.output}/intermediate",
    mode: params.publish_dir_mode,
    pattern: "*.gz",
    overwrite: false,
    saveAs: { fn ->
          fn.startsWith("keep_") ? "map-to-keep/${fn.replaceAll(~'^keep_', '')}" : "map-to-remove/${fn}"
    }
  )

  if ( !params.keep ) {
    publishDir (
      path: params.output,
      overwrite: false,
      mode: params.publish_dir_mode,
      pattern: "*.gz",
      saveAs: { fn ->
            fn.endsWith('.unmapped.fastq.gz') ? "clean/${fn}".replaceAll(~'.unmapped.fastq.gz$', '.fastq.gz') :
            fn.endsWith('.mapped.fastq.gz') ? "removed/${fn}".replaceAll(~'.mapped.fastq.gz$', '.fastq.gz') :
            fn.endsWith('.unmapped.fastq.gz') ? "clean/${fn}".replaceAll(~'.unmapped.fastq.gz$', '.fastq.gz') :
            fn.endsWith('.mapped.fastq.gz') ? "removed/${fn}".replaceAll(~'.mapped.fastq.gz$', '.fastq.gz') :
            fn
      }
    )
  }

  input:
  tuple val(name), val(type), path(bam)

  output:
  tuple val(name), val(type), path('*.fast*.gz')

  script:
  if ( params.lib_pairedness == 'paired' ) {
    """
    samtools fastq -@ ${task.cpus} -1 ${bam.baseName}_1.fastq -2 ${bam.baseName}_2.fastq -s ${bam.baseName}_singleton.fastq ${bam}
    gzip --no-name *.fastq
    """
  } else if ( params.lib_pairedness == 'single' ) {
    dtype = (params.input_type == 'fasta') ? 'a' : 'q'
    """
    samtools fastq -@ ${task.cpus} -0 ${bam.baseName}.fast${dtype} ${bam}
    gzip --no-name *.fast${dtype}
    """
  } else {
    error "Invalid pairedness: ${params.lib_pairedness}"
  }
  stub:
  dtype = (params.input_type == 'fasta') ? 'a' : 'q'
  """
  touch ${bam.baseName}_1.fast${dtype}.gz ${bam.baseName}_2.fast${dtype}.gz
  """
}

process idxstats_from_bam {
  label 'minimap2'

  publishDir "${params.output}/qc", mode: 'copy', pattern: "${bam.baseName}.idxstats.tsv", overwrite: false, enabled: false

  input:
  tuple val(name), val(type), path(bam), path(bai)

  output:
  tuple val(name), val(type), path('*.idxstats.tsv')

  script:
  """
  samtools idxstats ${bam} > ${bam.baseName}.idxstats.tsv
  """
  stub:
  """
  touch ${bam.baseName}.idxstats.tsv
  """
}

process flagstats_from_bam {
  label 'minimap2'

  publishDir "${params.output}/qc", mode: params.publish_dir_mode, pattern: "${bam.baseName}.flagstats.txt", overwrite: false, enabled: false

  input:
  tuple val(name), val(type), path(bam), path(bai)

  output:
  tuple val(name), val(type), path('*.flagstats.txt')

  script:
  """
  samtools flagstats ${bam} > ${bam.baseName}.flagstats.txt
  """
  stub:
  """
  touch ${bam.baseName}.flagstats.txt
  """
}

process sort_bam {
  label 'minimap2'

  input:
  tuple val(name), val(type), path(bam)

  output:
  tuple val(name), val(type), path("${bam.baseName}.sorted.bam")

  script:
  """
  mv ${bam} ${bam}.tmp
  samtools sort -@ ${task.cpus} ${bam}.tmp > ${bam.baseName}.sorted.bam
  """
  stub:
  """
  touch ${bam.baseName}.bam
  """
}

process index_bam {
  label 'minimap2'

  publishDir (
    path: "${params.output}/intermediate",
    mode: params.publish_dir_mode,
    pattern: "*.sorted.bam{,.bai}",
    overwrite: false,
    saveAs: { fn ->
          fn.startsWith("keep_") ? "map-to-keep/${fn.replaceAll(~'^keep_', '')}" : "map-to-remove/${fn}"
    }
  )

  input:
  tuple val(name), val(type), path(bam)

  output:
  tuple val(name), val(type), path(bam), path('*.bai')

  script:
  """
  samtools index -@ ${task.cpus} ${bam}
  """
  stub:
  """
  touch ${bam}.bai
  """
}

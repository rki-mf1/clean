process compress_reads {
  label 'basics'

  publishDir "${params.output}/${name}/${params.tool}", mode: 'copy', pattern: "*.gz"

  input:
  tuple val(name), val(type), path(reads)

  output:
  tuple val(name), val(type), path("*.fast{q,a}.gz")

  script:
  if ( params.mode == 'paired' ) {
    """
    pigz -fc -p ${task.cpus} ${reads[0]} > ${name}_1.${type}.fastq.gz 
    pigz -fc -p ${task.cpus} ${reads[1]} > ${name}_2.${type}.fastq.gz

    if [ -f "${reads[2]}" ]; then
      pigz -fc -p ${task.cpus} ${reads[2]} > ${name}.${type}_singleton.fastq.gz
    fi
    """
  } else if ( params.mode == 'single' ) {
    dtype = (params.seq_type == 'fasta') ? 'a' : 'q'
    """
    pigz -fc -p ${task.cpus} ${reads} > ${name}.${type}.fast${dtype}.gz
    """
  } else {
    error "Invalid mode: ${params.mode}"
  }
  stub:
  """
  touch ${name}.${type}.fasta.gz ${name}.${type}.fastq.gz
  """
}

process get_number_of_records {
  label 'smallTask'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), env(TOTALRECORDS), emit: TOTALRECORDS

  script:
  if ( params.mode == 'paired' ) {
    """
    if [[ ${reads[0]} =~ \\.gz\$ ]]; then
      TOTALRECORDS_1=\$(zcat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALRECORDS_2=\$(zcat ${reads[1]} | echo \$((`wc -l`/4)))
    else
      TOTALRECORDS_1=\$(cat ${reads[0]} | echo \$((`wc -l`/4)))
      TOTALRECORDS_2=\$(cat ${reads[1]} | echo \$((`wc -l`/4)))
    fi
    TOTALRECORDS=\$(( TOTALRECORDS_1+TOTALRECORDS_2 ))
    """
  } else if ( params.mode == 'single' && params.seq_type != 'fasta' ) {
    """
    if [[ ${reads} =~ \\.gz\$ ]]; then
      TOTALRECORDS=\$(zcat ${reads} | echo \$((`wc -l`/4)))
    else
      TOTALRECORDS=\$(cat ${reads} | echo \$((`wc -l`/4)))
    fi
    """
  } else if ( params.seq_type == 'fasta' ) {
    """
    if [[ ${reads} =~ \\.gz\$ ]]; then
      TOTALCONTIGS=\$(zgrep '^>' ${reads} | wc -l)
    else
      TOTALCONTIGS=\$(grep '^>' ${reads} | wc -l)
    fi
    """
  } else {
    error "Invalid mode: ${params.mode} or seq_type: ${params.seq_type}"
  }
  stub:
  """
  TOTALRECORDS=42
  """
}

process minimap2Stats {
  label 'smallTask'
  
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "stats.txt" 

  input:
  tuple val(name), path(idxstats), val (totalreads), val(ambiguousreads)

  output:
  tuple val(name), path ('stats.txt')
  path("${name}_minimap2_stats.tsv"), emit: tsv

  script:

  if ("${ambiguousreads}" == 'NULL'){
    header = 'Sample Name\tClean reads\tUnambiguous mapped reads'
    visible_val = ""
    val = '0'
  } else {
    header = 'Sample Name\tClean reads\tUnambiguous mapped reads\tAmbiguous mapped reads'
    visible_val = "\t${ambiguousreads}"
    val = "${ambiguousreads}"
  }

  """
  MAPPEDSUM=\$(awk -F '\\t' '{sum += \$3} END {print sum}' idxstats.tsv)
  UNPROMAPPEDSUM=\$(awk -F '\\t' '/^[^*]/ {sum += \$4} END {print sum}' idxstats.tsv)
  PROPMAP=\$((\$MAPPEDSUM-\$UNPROMAPPEDSUM))

  MAP=\$(awk -v map=\$MAPPEDSUM -v unpromap=\$UNPROMAPPEDSUM -v tot=${totalreads} 'BEGIN {perc=(map-unpromap)/tot*100; print map-unpromap " ("  perc " %) reads were properly mapped; of these:"}')

  FA=\$(awk -v tot=${totalreads} -F '\\t' '/^[^*]/ {propmap=\$3-\$4; if (propmap != 0) print "\\t\\t" propmap " (" propmap/tot*100  "%) reads aligned to " \$1}' idxstats.tsv;)

  touch stats.txt
  cat <<EOF >> stats.txt
  ${totalreads} reads in total; of these:
  \t\$MAP
  \$FA
  EOF

  touch ${name}_minimap2_stats.tsv
  cat <<EOF >> ${name}_minimap2_stats.tsv
  ${header}
  ${name}\t\$((${totalreads}-\$PROPMAP))\t\$((\$PROPMAP-${val}))${visible_val}
  EOF
  """
  stub:
  """
  touch stats.txt ${name}_minimap2_stats.tsv
  """
}

process bbdukStats {
  label 'smallTask'

  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "stats.txt"

  input:
  tuple val(name), path (bbdukStats)

  output:
  tuple val(name), path ('stats.txt')
  path("${name}_bbduk_stats.tsv"), emit: tsv

  script:
  """
  TOTAL=\$(grep '#Total' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MNUM=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MPER=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$3}')

  FA=\$(awk -F '\\t' '/^[^#]/ {print "\\t\\t"\$2" ("\$3") aligned to "\$1}' ${bbdukStats})

  touch stats.txt
  cat <<EOF >> stats.txt
  \$TOTAL reads in total; of these:
  \t\$MNUM (\$MPER) reads were properly mapped; of these:
  \$FA
  EOF

  touch ${name}_bbduk_stats.tsv
  cat <<EOF >> ${name}_bbduk_stats.tsv
  Sample Name\tClean reads\tMapped reads
  ${name}\t\$((\$TOTAL-\$MNUM))\t\$MNUM
  EOF
  """
  stub:
  """
  touch stats.txt ${name}_bbduk_stats.tsv
  """
}

process writeLog {
  label 'smallTask'

  publishDir "${params.output}/${name}/${params.tool}", mode: 'copy', pattern: "log.txt"
  
  input:
    tuple val(name), path (reads)
    val db

  output:
    path 'log.txt'
  
  script:
  """
  touch log.txt
  cat <<EOF >> log.txt
  Input reads:\t${reads}
  Contamination:\t${db}
  
  Statistics summary:\t${params.output}/${name}/${params.tool}/stats.txt
  EOF
  """
  stub:
  """
  touch log.txt
  """
}

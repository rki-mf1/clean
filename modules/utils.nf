process compress_reads {
  label 'basics'

  publishDir "${params.output}/${params.tool}", mode: 'copy', pattern: "*.gz"

  input:
  tuple val(name), val(type), path(reads)

  output:
  tuple val(name), val(type), path("*.fast{q,a}.gz")

  script:
  if ( params.lib_pairedness == 'paired' ) {
    """
    pigz -fc -p ${task.cpus} ${reads[0]} > ${name}_1.${type}.fastq.gz 
    pigz -fc -p ${task.cpus} ${reads[1]} > ${name}_2.${type}.fastq.gz

    if [ -f "${reads[2]}" ]; then
      pigz -fc -p ${task.cpus} ${reads[2]} > ${name}.${type}_singleton.fastq.gz
    fi
    """
  } else if ( params.lib_pairedness == 'single' ) {
    dtype = (params.input_type == 'fasta') ? 'a' : 'q'
    """
    pigz -fc -p ${task.cpus} ${reads} > ${name}.${type}.fast${dtype}.gz
    """
  } else {
    error "Invalid mode: ${params.lib_pairedness}"
  }
  stub:
  dtype = (params.input_type == 'fasta') ? 'a' : 'q'
  """
  touch ${name}.${type}.fast${dtype}.gz
  """
}

process get_number_of_records {
  label 'smallTask'

  input:
  tuple val(name), path(reads)

  output:
  tuple val(name), env(TOTALRECORDS), emit: TOTALRECORDS

  script:
  if ( params.lib_pairedness == 'paired' ) {
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
  } else if ( params.lib_pairedness == 'single' && params.input_type != 'fasta' ) {
    """
    if [[ ${reads} =~ \\.gz\$ ]]; then
      TOTALRECORDS=\$(zcat ${reads} | echo \$((`wc -l`/4)))
    else
      TOTALRECORDS=\$(cat ${reads} | echo \$((`wc -l`/4)))
    fi
    """
  } else if ( params.input_type == 'fasta' ) {
    """
    if [[ ${reads} =~ \\.gz\$ ]]; then
      TOTALCONTIGS=\$(zgrep '^>' ${reads} | wc -l)
    else
      TOTALCONTIGS=\$(grep '^>' ${reads} | wc -l)
    fi
    """
  } else {
    error "Invalid pairedness: ${params.lib_pairedness} or input_type: ${params.input_type}"
  }
  stub:
  """
  TOTALRECORDS=42
  """
}

process get_read_names {
  label 'minimap2'

  input:
  tuple val(name), val(type), path(bam)
  
  output:
  path("${name}_read_names.csv")
  
  script:
  """
  samtools view ${bam} | cut -f1 | sort | uniq > ${name}_read_names.csv
  """
}

process split_bam {
  label 'pysam'
  echo true

  input:
  path(read_name_list)
  tuple val(name), val(type), path(mapped_bam)

  output:
  tuple val(name), val(type), path('mapped.fq'), emit: mapped
  tuple val(name), val(type), path('unmapped.fq'), emit: unmapped

  script:
  """
  #!/usr/bin/env python3

  import pysam

  reads = set()
  with open('${read_name_list}', 'r') as infile:
    for line in infile:
      reads.add(line.strip())
  print(reads)
  
  # split bam into mapped (not in list)
  # and "unmapped" (in list)

  bamfile = pysam.AlignmentFile('nanopore.mapped.bam', 'rb')
  for read in bamfile.fetch(until_eof=True):
    if ( read.query_name in reads):
      #write to unmapped
    else:
      # write to mapped
  """
}

process bbdukStats {
  label 'smallTask'

  publishDir "${params.output}/bbduk", mode: 'copy', pattern: "${name}_stats.txt"

  input:
  tuple val(name), path (bbdukStats)

  output:
  tuple val(name), path ("${name}_stats.txt")
  path("${name}_bbduk_stats.tsv"), emit: tsv

  script:
  """
  TOTAL=\$(grep '#Total' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MNUM=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MPER=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$3}')

  FA=\$(awk -F '\\t' '/^[^#]/ {print "\\t\\t"\$2" ("\$3") aligned to "\$1}' ${bbdukStats})

  touch ${name}_stats.txt
  cat <<EOF >> ${name}_stats.txt
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
  touch ${name}_stats.txt ${name}_bbduk_stats.tsv
  """
}

process writeLog {
  label 'smallTask'

  publishDir "${params.output}/${params.tool}", mode: 'copy', pattern: "log.txt"
  
  input:
    val db
    path (reads)

  output:
    path 'log.txt'
  
  script:

  """
  touch log.txt
  cat <<EOF >> log.txt
  Input reads:\t${reads}
  Contamination:\t${db}
  EOF
  """
  stub:
  """
  touch log.txt
  """
}

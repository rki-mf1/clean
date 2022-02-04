process fastqc {
  label 'fastqc'

  input:
  tuple val(name), val(type), path(reads)

  output:
  tuple val(name), val(type), path("*_fastqc.zip"), emit: zip

  script:
  """
  fastqc --noextract -t ${task.cpus} ${reads}
  """
  stub:
  """
  touch ${name}_fastqc.zip
  """
}

process nanoplot {
  label 'nanoplot'
  errorStrategy { task.exitStatus in 1 ? 'ignore' : 'terminate' }

  input:
  tuple val(name), val(type), path(reads)

  output:
  tuple val(name), path("*.html"), path("*.png"), path("*.pdf"), file("${name}_${type}_read_quality.txt")
  tuple val(name), val(type), path("${name}_${type}_read_quality_report.html"), emit: html
  
  script:
  """
  NanoPlot -t ${task.cpus} --fastq ${reads} --title ${name}_${type} --color darkslategrey --N50 --plots hex --loglength -f png --store
  NanoPlot -t ${task.cpus} --pickle NanoPlot-data.pickle --title ${name}_${type} --color darkslategrey --N50 --plots hex --loglength -f pdf
  mv NanoPlot-report.html ${name}_${type}_read_quality_report.html
  mv NanoStats.txt ${name}_${type}_read_quality.txt
  """
  stub:
  """
  touch ${name}_${type}_read_quality_report.html ${name}_${type}_read_quality.txt fuu.png fuu.pdf
  """
}

process format_nanoplot_report {
    label 'smallTask'
    
    input:
    tuple val(name), val(type), path(nanoplot_report)

    output:
    path("*_mqc.html")

    script:
    """
    sed -e '25,30d;34,45d' ${nanoplot_report} > ${nanoplot_report}.tmp
    echo "<!--" > tmp
    echo "id: 'nanoplot_${name}_${type}'" >> tmp
    echo "section_name: 'NanoPlot: ${name}, ${type}'" >> tmp
    echo "-->"  >> tmp
    cat tmp ${nanoplot_report}.tmp > ${nanoplot_report.baseName}_mqc.html
    rm -f *tmp
    """
    stub:
    """
    touch ${nanoplot_report.baseName}_mqc.html
    """
}

process quast {
  label 'quast'
  errorStrategy { task.exitStatus in 4 ? 'ignore' : 'terminate' }

  input:
  tuple val(name), val(type), path(fasta)

  output:
  path("${name}_${type}_report.tsv"), emit: report_tsv
  path("quast_${name}_${type}")

  script:
  """
  quast.py -o quast_${name}_${type} -t ${task.cpus} ${fasta}
  cp quast_${name}_${type}/report.tsv ${name}_${type}_report.tsv
  """
  stub:
  """
  touch ${name}_${type}_report.tsv quast_${name}_${type}
  """
}

process multiqc {
  label 'multiqc'
  label 'smallTask'
  
  publishDir "${params.output}/${params.multiqc_dir}", pattern: 'multiqc_report.html'
  
  input:
  path(config)
  path(fastqc)
  path(nanoplot)
  path(quast)
  path(mapping_stats)
  path(idxstats)
    
  output:
  path "multiqc_report.html"
  
  script:
  """
  multiqc . -s -c ${config}
  """
  stub:
  """
  touch multiqc_report.html
  """
}

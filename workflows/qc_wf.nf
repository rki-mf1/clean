include { fastqc; nanoplot; format_nanoplot_report; quast; multiqc } from '../modules/qc'

workflow qc {
  take:
    input
    input_type
    bbduk_summary
    idxstats
    flagstats
    multiqc_config
  main:
    if ( input_type == 'fasta' ){
      quast(input)
      report = quast.out.report_tsv
    } else if ( input_type == 'nano' ) {
      nanoplot(input)
      format_nanoplot_report(nanoplot.out.html)
      report = format_nanoplot_report.out
    } else if ( input_type.contains('illumina') ){
      fastqc(input)
      report = fastqc.out.zip.map{ it -> it[-1] }
    } else { error "Invalid input type: ${input_type}" }

    multiqc(multiqc_config, report.collect(), bbduk_summary.collect().ifEmpty([]), idxstats.map{ it -> it[1] }.collect().ifEmpty([]), flagstats.map{ it -> it[1] }.collect().ifEmpty([]))
}

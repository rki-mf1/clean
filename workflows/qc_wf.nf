include { fastqc; nanoplot; format_nanoplot_report; quast; multiqc } from '../modules/qc'

workflow qc_fasta {
  take:
    fasta_input
    fasta_output
  main:
    quast(fasta_input.concat(fasta_output))
  emit:
    quast.out.report_tsv
}

workflow qc_nano {
  take:
    nano_input
    nano_output
  main:
    nanoplot(nano_input.concat(nano_output))
    format_nanoplot_report(nanoplot.out.html)
  emit:
    format_nanoplot_report.out
}

workflow qc_illumina {
  take:
    illumina_input
    illumina_output
  main:
    fastqc(illumina_input.concat(illumina_output))
  emit:
    fastqc.out.zip.map{ it -> it[-1] }
}

workflow qc{
  take:
    multiqc_config
    fastqc
    nanoplot
    quast
    mapping_stats
    idxstats
  main:
    multiqc(multiqc_config, fastqc, nanoplot, quast, mapping_stats, idxstats)
}
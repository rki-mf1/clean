include { minimap2 } from '../modules/minimap2'
include { bbduk } from '../modules/bbmap'
include { compress_reads; writeLog ; bbdukStats } from '../modules/utils'
include { split_bam; fastq_from_bam ; idxstats_from_bam ; flagstats_from_bam ; index_bam ; sort_bam } from '../modules/alignment_processing'

workflow clean {
    take:
        input
        contamination

    main:
        if ( params.bbduk ) {
            // map
            bbduk(input, contamination)
            // compress reads
            compress_reads(bbduk.out.cleaned_reads.concat(bbduk.out.contaminated_reads))
            // log & stats
            writeLog(contamination, input.map{ it -> it[1] }.collect())
            bbdukStats(bbduk.out.stats)
            // define output
            bbduk_summary = bbdukStats.out.tsv
            idxstats = []
            flagstats = []
            out_reads = bbduk.out.cleaned_reads.concat(bbduk.out.contaminated_reads)
        } 
        else {
            minimap2(input, contamination) | sort_bam | index_bam | ( idxstats_from_bam & flagstats_from_bam )

            split_bam(minimap2.out.bam)
            contamination_bam = split_bam.out.mapped
            cleaned_bam = split_bam.out.unmapped
        //     if ( params.control && 'dcs' in params.control.split(',') ) {
        //         filter_true_dcs_alignments(contamination_bam)
        //         contamination_bam = concat_bams(filter_true_dcs_alignments.out.no_dcs.mix(filter_true_dcs_alignments.out.true_dcs))
        //         cleaned_bam = concat_bams(cleaned_bam.mix(filter_true_dcs_alignments.out.false_dcs))
        //     }
        //     if ( params.min_clip ) {
        //         filter_soft_clipped_alignments(contamination_bam, params.min_clip)
        //         contamination_bam = filter_soft_clipped_alignments.out.bam_unam
        //         cleaned_bam = concat_bams(cleaned_bam.mix(filter_soft_clipped_alignments.out.bam_am)
        //     }
            fastq_from_bam(contamination_bam.mix(cleaned_bam))
            // compress reads
            compress_reads(fastq_from_bam.out)
            // log & stats
            writeLog(contamination, input.map{ it -> it[1] }.collect())
            // define output
            bbduk_summary = []
            idxstats = idxstats_from_bam.out
            flagstats = flagstats_from_bam.out
            out_reads = compress_reads.out
        }

    emit:
        bbduk_summary
        idxstats
        flagstats
        out_reads
}
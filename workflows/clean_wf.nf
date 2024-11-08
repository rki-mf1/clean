include { minimap2 } from '../modules/minimap2'
include { bwa_index; bwa } from '../modules/bwa'
include { bbduk } from '../modules/bbmap'
include { bbduk_stats } from '../modules/utils'
include { split_bam; fastq_from_bam ; idxstats_from_bam ; flagstats_from_bam ; index_bam as index_bam; index_bam as index_bam2; sort_bam ; filter_true_dcs_alignments ; merge_bam as merge_bam1 ; merge_bam as merge_bam2 ; merge_bam as merge_bam3 ; merge_bam as merge_bam4 ; filter_soft_clipped_alignments } from '../modules/alignment_processing'

workflow clean {
    take:
        input
        contamination
        dcs_ends_bed
        map_target

    main:
        if ( params.bbduk ) {
            // map
            bbduk(input, contamination, map_target)
            bbduk_stats(bbduk.out.stats)
            // define output
            bbduk_summary = bbduk_stats.out.tsv
            idxstats = Channel.empty()
            flagstats = Channel.empty()
            out_reads = bbduk.out.cleaned_reads.concat(bbduk.out.contaminated_reads)
            bams_bai = Channel.empty()
            sort_bam_ch = Channel.empty()
        }
        else {
            if ( params.bwa ) {
                bwa_index(contamination)
                bwa(input, bwa_index.out, contamination) | sort_bam | index_bam | ( idxstats_from_bam & flagstats_from_bam )
                split_bam(bwa.out.bam)
            } else {
                minimap2(input, contamination) | sort_bam | index_bam | ( idxstats_from_bam & flagstats_from_bam )
                split_bam(minimap2.out.bam)
            }
            contamination_bam = split_bam.out.mapped
            cleaned_bam = split_bam.out.unmapped
            if ( params.control && 'dcs' in params.control.split(',') && params.dcs_strict ) {
                filter_true_dcs_alignments(contamination_bam.map{ it -> [it[0], it[2]] }, dcs_ends_bed)
                
                contamination_bam = merge_bam1(contamination_bam.mix(filter_true_dcs_alignments.out.true_dcs).mix(filter_true_dcs_alignments.out.no_dcs).groupTuple(by: [0,1]))
                cleaned_bam = merge_bam2(cleaned_bam.mix(filter_true_dcs_alignments.out.false_dcs).groupTuple(by: [0,1]))
            }
            if ( params.min_clip ) {
                filter_soft_clipped_alignments(contamination_bam.map{ it -> tuple(it[0], it[2]) }, params.min_clip)

                contamination_bam = merge_bam3(contamination_bam.mix(filter_soft_clipped_alignments.out.bam_ok_clipped).groupTuple(by: [0,1]))
                cleaned_bam = merge_bam4(cleaned_bam.mix(filter_soft_clipped_alignments.out.bam_clipped).groupTuple(by: [0,1]))
            }
            index_bam2(contamination_bam.mix(cleaned_bam))
            bams_bai = index_bam2.out

            fastq_from_bam(contamination_bam.mix(cleaned_bam))

            // define output
            bbduk_summary = Channel.empty()
            idxstats = idxstats_from_bam.out
            flagstats = flagstats_from_bam.out
            out_reads = fastq_from_bam.out
	          sort_bam_ch = sort_bam.out
        }

    emit:
        bbduk_summary
        idxstats
        flagstats
        out_reads
        bams_bai
	      sort_bam_ch
}

include { minimap2 } from '../modules/minimap2' addParams( mode: lib_type; seq_type: seq_type )
include { bbduk } from '../modules/bbmap' addParams( mode: lib_type )

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
            writeLog(input, contamination)
            bbdukStats(bbduk.out.stats)
            stats = bbdukStats.out.tsv  
        } else {
            minimap2(input_ch, contamination)
            split_sam(minimap2.out.sam)
            contamination_bam = split_sam.out.mapped
            cleaned_bam = split_sam.out.unmapped
            if ( params.control && 'dcs' in params.control.split(',') ) {
                filter_true_dcs_alignments(contamination_bam)
                contamination_bam = concat_bams(filter_true_dcs_alignments.out.no_dcs.mix(filter_true_dcs_alignments.out.true_dcs))
                cleaned_bam = concat_bams(cleaned_bam.mix(filter_true_dcs_alignments.out.false_dcs))
            }
            if ( params.min_clip ) {
                filter_soft_clipped_alignments(contamination_bam, params.min_clip)
                contamination_bam = filter_soft_clipped_alignments.out.bam_unam
                cleaned_bam = concat_bams(cleaned_bam.mix(filter_soft_clipped_alignments.out.bam_am)
            }
            fastq_from_bam(contamination_bam)
            fastq_from_bam(cleaned_bam)
            // log & stats
            writeLog(input, contamination)
            get_number_of_records(input)
            // minimap2Stats(make_mapped_bam.out.idxstats.join(get_number_of_records.out).join( number_ambiguous_reads_ch ) )
            // stats = minimap2Stats.out.tsv
        }

    // emit:
}
include { clean as keep_map } from './clean_wf'

include { get_read_names; filter_fastq_by_name } from '../modules/utils'
include { idxstats_from_bam ; flagstats_from_bam } from '../modules/alignment_processing'

workflow keep {
    take:
        input
        keep_reference
        dcs_ends_bed
        mapped_clean_fastq
        unmapped_clean_fastq

    main:
        keep_map(input.map{ it -> ['keep_'+it[0], it[1]]}, keep_reference, dcs_ends_bed)
        
        keep_reads_bam = keep_map.out.contamination_bam_bai
        get_read_names(keep_reads_bam.map{it -> [it[0], it[1]]})
        // works also for multiple samples?
        
        // mv keep mapped reads from cleaned mapped to clean unmapped
        filter_fastq_by_name(get_read_names.out, mapped_clean_fastq.join(unmapped_clean_fastq))

        //  QC
        idxstats = idxstats_from_bam(keep_reads_bam)
        flagstats = flagstats_from_bam(keep_reads_bam)

        idxstats = idxstats_from_bam.out
        flagstats = flagstats_from_bam.out
        out_reads = filter_fastq_by_name.out.mapped_no_keep.mix(filter_fastq_by_name.out.unmapped_keep)
    emit:
        idxstats
        flagstats
        out_reads
        keep_reads_bam
}
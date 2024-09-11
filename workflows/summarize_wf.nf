include { BAM_STATISTICS; COMBINE_BAM_STATISTICS } from '../modules/summarize'

workflow summarize {
    
    take:
       sorted_bam
           
    main:
       BAM_STATISTICS(sorted_bam)
       
       collected_bam_statistics_ch = BAM_STATISTICS.out.collect( sort: {a, b -> a[0].getBaseName() <=> b[0].getBaseName()} )
	
	
	COMBINE_BAM_STATISTICS(collected_bam_statistics_ch)
	
	// define output
        bamstats = BAM_STATISTICS.out
        summarystats = COMBINE_BAM_STATISTICS.out
        
    emit:
        bamstats
        summarystats
}

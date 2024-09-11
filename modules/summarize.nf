process BAM_STATISTICS {
   
    label 'minimap2'
    
    publishDir (
    
    	path: "${params.output}/intermediate",
	mode: params.publish_dir_mode,
	pattern: "*.mappedstats.txt",
	enabled: !params.no_intermediate,
	saveAs: { fn ->
          fn.startsWith("keep_") ? "map-to-keep/${fn.replaceAll(~'^keep_', '')}" : "map-to-remove/${fn}"
    }
    )
    
    input:
    tuple val(name), val(type), path(bam)
    

    output:
    path("*.mappedstats.txt")

    
    script:
    """ 
    bam_stats.sh ${name}

    """
    
    stub:
    """
    touch ${bam.baseName}.mappedstats.txt
    """
}


process COMBINE_BAM_STATISTICS {
    
    label 'minimap2'
    
    publishDir (
    
    	path: "${params.output}/intermediate",
	mode: params.publish_dir_mode,
	pattern: "combined_bam_stats.txt",
	enabled: !params.no_intermediate,
	saveAs: { fn ->
          fn.startsWith("keep_") ? "map-to-keep/${fn.replaceAll(~'^keep_', '')}" : "map-to-remove/${fn}"
    }
    )
    
    tag { 'combine bam statistics files'} 
    
    
    input:
    path(bam_statistics_files)
    

    output:
    path("combined_bam_stats.txt")

    
    script:
    """ 
    BAM_STATISTICS_FILES=(${bam_statistics_files})
    
    for index in \${!BAM_STATISTICS_FILES[@]}; do
    BAM_STATISTICS_FILE=\${BAM_STATISTICS_FILES[\$index]}
    
    # add header line if first file
    if [[ \$index -eq 0 ]]; then
      echo "\$(head -1 \${BAM_STATISTICS_FILE})" >> combined_bam_stats.txt
    fi
    echo "\$(awk 'FNR==2 {print}' \${BAM_STATISTICS_FILE})" >> combined_bam_stats.txt
    done

    """
    
    stub:
    """
    touch combined_bam_stats.txt
    """
}

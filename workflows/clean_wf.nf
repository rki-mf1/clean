include { concat_contamination } from '../modules/get_host' addParams( tool: tool )

include { minimap2 } from '../modules/minimap2' addParams( mode: lib_type; seq_type: seq_type )
include { bbduk } from '../modules/bbmap' addParams( mode: lib_type )

workflow clean {
    take:
        input
        contamination

    main:
        concat_contamination( input.map{ it -> it[0] }, contamination )
        if ( params.bbduk ) {
            // map
            bbduk(input, concat_contamination.out.fa)
        } else {
            minimap2(input_ch, concat_contamination.out.fa)
        }

    // emit:
}
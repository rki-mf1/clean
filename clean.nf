#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Nextflow -- Decontamination Pipeline
Author: marie.lataretu@uni-jena.de
Author: hoelzer.martin@gmail.com
*/

// Parameters sanity checking

Set valid_params = ['max_cores', 'cores', 'max_memory', 'memory', 'profile', 'help', 'input', 'input_type', 'list', 'host', 'own', 'control', 'keep', 'rm_rrna', 'bwa', 'bbduk', 'bbduk_kmer', 'bbduk_qin', 'reads_rna', 'min_clip', 'dcs_strict', 'output', 'multiqc_dir', 'nf_runinfo_dir', 'databases', 'cleanup_work_dir','condaCacheDir', 'singularityCacheDir', 'singularityCacheDir', 'cloudProcess', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process', 'publish_dir_mode', 'no_intermediate', 'skip_qc'] // don't ask me why there is also 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'
def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
}
if (params.input.contains('.clean.') ) {
  exit 1, "ERROR: Input files cannot contain `.clean.`\n"
}

/**************************
* META & HELP MESSAGES
**************************/

/*
Comment section: First part is a terminal print for additional user information,
followed by some help statements (e.g. missing input) Second part is file
channel input. This allows via --list to alter the input of --nano & --illumina
to add csv instead. name,path or name,pathR1,pathR2 in case of illumina
*/

// terminal prints
if (params.help) { exit 0, helpMSG() }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $nextflow.timestamp"
println "Output directory name:"
println "  $params.output"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
println "Database location:"
println "  $params.databases"
if ( workflow.profile.contains('singularity') ) {
    println "Singularity cache directory:"
    println "  $params.singularityCacheDir"
}
if ( workflow.profile.contains('conda') ) {
    println "Conda cache directory:"
    println "  $params.condaCacheDir"
}
println "Configuration files:"
println "  $workflow.configFiles"
println "Cmd line:"
println "  $workflow.commandLine\u001B[0m"
if (workflow.repository != null){ println "\033[2mGit info: $workflow.repository - $workflow.revision [$workflow.commitId]\u001B[0m" }
println " "
if (workflow.profile == 'standard' || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println "\033[2mMemory to use: $params.memory, maximal memory to use: $params.max_memory\u001B[0m"
    println " "
}
if ( !workflow.revision ) {
    println "\033[0;33mWARNING: Not a stable execution. Please use -r for full reproducibility.\033[0m\n"
}
def folder = new File(params.output)
if ( folder.exists() ) {
    println "\033[0;33mWARNING: Output folder already exists. Results might be overwritten! You can adjust the output folder via [--output]\033[0m\n"
}
if ( workflow.profile.contains('singularity') ) {
    println "\033[0;33mWARNING: Singularity image building sometimes fails!"
    println "Multiple resumes (-resume) and --max_cores 1 --cores 1 for local execution might help.\033[0m\n"
}

Set controls = ['phix', 'dcs', 'eno']
Set hosts = ['hsa', 'mmu', 'cli', 'csa', 'gga', 'eco', 'sc2', 't2t']
Set input_types = ['nano', 'illumina', 'illumina_single_end', 'fasta']

if ( params.profile ) { exit 1, "--profile is wrong, use -profile" }
if ( params.input == '' || !params.input_type == '' ) { exit 1, "Missing required input parameters [--input] and [--input_type]" }

if ( params.input_type ) { if ( ! (params.input_type in input_types ) ) { exit 1, "Choose one of the the input types with --input_type: " + input_types } }

if ( params.control ) { for( String ctr : params.control.split(',') ) if ( ! (ctr in controls ) ) { exit 1, "Wrong control defined (" + ctr + "), use one of these: " + controls } }
if ( params.input_type == 'nano' && params.control && 'dcs' in params.control.split(',') && 'eno' in params.control.split(',') ) { exit 1, "Please choose either eno (for ONT dRNA-Seq) or dcs (for ONT DNA-Seq)." }
if ( params.host ) { for( String hst : params.host.split(',') ) if ( ! (hst in hosts ) ) { exit 1, "Wrong host defined (" + hst + "), use one of these: " + hosts } }
if ( !params.host && !params.own && !params.control && !params.rm_rrna ) { exit 1, "Please provide a control (--control), a host tag (--host), a FASTA file (--own) or set --rm_rrna for rRNA removal for the clean up."}

/**************************
* INPUT CHANNELS
**************************/

if ( params.input_type == 'illumina' ) {
  if ( params.list ) { input_ch = Channel
    .fromPath( params.input, checkIfExists: true )
    .splitCsv()
    .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  } else { input_ch = Channel
      .fromFilePairs( params.input , checkIfExists: true )
  }
} else {
  if ( params.list ) {
    input_ch = Channel
     .fromPath( params.input, checkIfExists: true )
     .splitCsv()
     .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  } else { input_ch = Channel
    .fromPath( params.input, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
  }
}

// load control fasta sequence
if ( params.control ) {
  if ( 'phix' in params.control.split(',') ) {
    illuminaControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/phix.fa.gz' , checkIfExists: true )
    nanoControlBedChannel = []
  } else { illuminaControlFastaChannel = Channel.empty() }
  if ( 'dcs' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/dcs.fa.gz' , checkIfExists: true )
    nanoControlBedChannel = Channel.fromPath( workflow.projectDir + '/data/controls/dcs_artificial_ends.bed' , checkIfExists: true )
  } else if ( 'eno' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/eno.fa.gz' , checkIfExists: true )
    nanoControlBedChannel = []
  } else { nanoControlFastaChannel = Channel.empty() }
} else {
  nanoControlFastaChannel = Channel.empty()
  illuminaControlFastaChannel = Channel.empty()
  nanoControlBedChannel = []
}

// load rRNA DB
if ( params.rm_rrna ){
  rRNAChannel = Channel.fromPath( workflow.projectDir + '/data/rRNA/*.fasta.gz', checkIfExists: true )
} else{
  rRNAChannel = Channel.empty()
}

if ( params.host ) {
  hostNameChannel = Channel.from( params.host ).splitCsv().flatten()
} else {
  hostNameChannel = Channel.empty()
}

// user defined fasta sequence
if ( params.own && params.list ) {
  ownFastaChannel = Channel
    .fromPath( params.own, checkIfExists: true)
    .splitCsv().flatten().map{ it -> file( it, checkIfExists: true ) }
} else if ( params.own ) {
  ownFastaChannel = Channel
    .fromPath( params.own, checkIfExists: true)
} else {
  ownFastaChannel = Channel.empty()
}

// user defined fasta sequence to keep
if ( params.keep && params.list ) {
  keepFastaChannel = Channel
    .fromPath( params.keep, checkIfExists: true)
    .splitCsv().flatten().map{ it -> file( it, checkIfExists: true ) }
} else if ( params.keep ) {
  keepFastaChannel = Channel
    .fromPath( params.keep, checkIfExists: true)
}

multiqc_config = Channel.fromPath( workflow.projectDir + '/assets/multiqc_config.yml', checkIfExists: true )

tool = params.bbduk ? 'bbduk' : 'minimap2'
lib_pairedness = params.input_type == 'illumina' ? 'paired' : 'single'

/**************************
* MODULES
**************************/

include { prepare_contamination } from './workflows/prepare_contamination_wf' addParams( tool: tool )
include { check_own as prepare_keep } from './modules/prepare_contamination'
include { clean } from './workflows/clean_wf' addParams( tool: tool, lib_pairedness: lib_pairedness )
include { keep } from './workflows/keep_wf' addParams( tool: tool, lib_pairedness: lib_pairedness )
include { summarize } from './workflows/summarize_wf'
include { qc } from './workflows/qc_wf'

/**************************
* WORKFLOW ENTRY POINT
**************************/

workflow {
  prepare_contamination(nanoControlFastaChannel, illuminaControlFastaChannel, rRNAChannel, hostNameChannel, ownFastaChannel)
  contamination = prepare_contamination.out

  clean(input_ch, contamination, nanoControlBedChannel, 'map-to-remove')
    
  if (!params.bbduk) {
    bam_sort = clean.out.sort_bam_ch
    summarize(bam_sort)
  }
 
  if (params.keep){
    prepare_keep(keepFastaChannel)
    keep_fasta = prepare_keep.out

    mapped = clean.out.out_reads.filter{ it[1] == 'mapped' }
    unmapped = clean.out.out_reads.filter{ it[1] == 'unmapped' }

    un_mapped_clean_fastq = mapped.join(unmapped)

    keep(input_ch.map{ it -> ['keep_'+it[0], it[1]]}, keep_fasta.collect(), nanoControlBedChannel, un_mapped_clean_fastq)

  }

  if (!params.skip_qc) {
    qc(input_ch.map{ it -> tuple(it[0], 'input', it[1]) }.mix(clean.out.out_reads), params.input_type, clean.out.bbduk_summary, clean.out.idxstats, clean.out.flagstats, multiqc_config)
  }
}

/**************************
* --help
**************************/
def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ____________________________________________________________________________________________

    Workflow: Decontamination

    Clean your Illumina, Nanopore or any FASTA-formated sequence date. The output are the clean
    and as contaminated identified sequences. Per default minimap2 is used for aligning your sequences
    to a host but we recommend using the ${c_dim}--bbduk${c_reset} flag to switch to bbduk to clean short-read data.

    Use the ${c_dim}--host${c_reset} and ${c_dim}--control${c_reset} flag to download a host database or specify your ${c_dim}--own${c_reset} FASTA.

    ${c_yellow}Usage example:${c_reset}
    nextflow run rki-mf1/clean --input_type nano --input '*/*.fastq' --host eco --control dcs
    or
    nextflow run rki-mf1/clean --input_type illumina --input '*/*.R{1,2}.fastq' --own some_host.fasta --bbduk
    or
    nextflow run rki-mf1/clean --input_type illumina --input 'test/illumina*.R{1,2}.fastq.gz' --nano data/nanopore.fastq.gz --fasta data/assembly.fasta --host eco --control phix

    ${c_yellow}Input:${c_reset}
    ${c_green}--input_type nano                --input${c_reset} '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green}--input_type illumina            --input${c_reset} '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green}--input_type illumina_single_end --input${c_reset} '*.fastq.gz'                -> one sample per file
    ${c_green}--input_type fasta               --input${c_reset} '*.fasta.gz'                -> one sample per file
    ${c_dim} ...read above input from csv files:${c_reset} ${c_green}--list ${c_reset}
                         ${c_dim}required format: name,path for --input_type nano and --input_type fasta; name,pathR1,pathR2 for --illumina input_type; name,path for --input_type illumina_single_end${c_reset}

    ${c_yellow}Decontamination options:${c_reset}
    ${c_green}--host${c_reset}         Comma separated list of reference genomes for decontamination, downloaded based on this parameter [default: $params.host]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly, incl. mtDNA]
                                        - t2t [T2T Consortium: human genome w/ additional 200 Mbp, closed gaps, and more complete Y (T2T-CHM13+Yv2.0), incl. mtDNA]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly, incl. mtDNA]
                                        - csa [NCBI: GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic, incl. mtDNA]
                                        - gga [NCBI: Gallus_gallus.GRCg6a.dna.toplevel, incl. mtDNA]
                                        - cli [NCBI: GCF_000337935.1_Cliv_1.0_genomic, incl. mtDNA]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]
                                        - sc2 [ENA: MN908947.3 (Wuhan-Hu-1 complete genome)]${c_reset}
    ${c_green}--control${c_reset}       Comma separated list of common controls used in Illumina or Nanopore sequencing [default: $params.control]
                                        ${c_dim}Currently supported are:
                                        - phix [Illumina: enterobacteria_phage_phix174_sensu_lato_uid14015, NC_001422]
                                        - dcs [ONT DNA-Seq: a positive control (3.6 kb standard amplicon mapping the 3' end of the Lambda genome)]
                                        - eno [ONT RNA-Seq: a positive control (yeast ENO2 Enolase II of strain S288C, YHR174W)]${c_reset}
    ${c_green}--own ${c_reset}          Use your own FASTA sequences (comma separated list of files) for decontamination, e.g. host.fasta.gz,spike.fasta [default: $params.own]
    ${c_green}--keep ${c_reset}         Use your own FASTA sequences (comma separated list of files) to explicitly keep mapped reads, e.g. target.fasta.gz,important.fasta [default: $params.keep]
                                        Reads are assigned to a combined index for decontamination and keeping. The use of this parameter can prevent
                                        false positive hits and the accidental removal of reads due to (poor quality) mappings.
    ${c_green}--rm_rrna ${c_reset}      Clean your data from rRNA [default: $params.rm_rrna]
    ${c_green}--bwa${c_reset}           Add this flag to use BAW MEM instead of minimap2 for decontamination of short reads [default: $params.bwa]
    ${c_green}--bbduk${c_reset}         Add this flag to use bbduk instead of minimap2 for decontamination of short reads [default: $params.bbduk]
    ${c_green}--bbduk_kmer${c_reset}    Set kmer for bbduk [default: $params.bbduk_kmer]
    ${c_green}--bbduk_qin${c_reset}     Set quality ASCII encoding for bbduk [default: $params.bbduk_qin; options are: 64, 33, auto]
    ${c_green}--reads_rna${c_reset}     Add this flag for noisy direct RNA-Seq Nanopore data [default: $params.reads_rna]

    ${c_green}--min_clip${c_reset}      Filter mapped reads by soft-clipped length (left + right). If >= 1 total number; if < 1 relative to read length
    ${c_green}--dcs_strict${c_reset}    Filter out alignments that cover artificial ends of the ONT DCS to discriminate between Lambda Phage and DCS
    ${c_green}--skip_qc${c_reset}       Skip quality control steps (fastqc, nanoplot, multiqc, etc.) [default: $params.skip_qc]

    ${c_yellow}Compute options:${c_reset}
    --cores             Max cores per process for local use [default $params.cores]
    --max_cores         Max cores used on the machine for local use [default $params.max_cores]
    --memory            Max memory for local use, enter in this format '8.GB' [default: $params.memory]
    --output            Name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    CPU / RAM usage (may cause errors)
    -with-dag chart.html     Generates a flowchart for the process tree
    -with-timeline time.html Timeline (may cause errors)

    ${c_yellow}Computing:${c_reset}
    In particular for execution of the workflow on a HPC (LSF, SLURM) adjust the following parameters:
    --databases             Defines the path where databases are stored [default: $params.databases]
    --condaCacheDir         Defines the path where environments (conda) are cached [default: $params.condaCacheDir]
    --singularityCacheDir   Defines the path where images (singularity) are cached [default: $params.singularityCacheDir]

    ${c_yellow}Miscellaneous:${c_reset}
    --cleanup_work_dir      Deletes all files in the work directory after a successful completion of a run [default: $params.cleanup_work_dir]
                            ${c_dim}warning: if true, the option will prevent the use of the resume feature!${c_reset}
    --no_intermediate       Do not save intermediate .bam/fastq/etc files into the `results/intermediate/` directory [default: $params.cleanup_work_dir]
                            Saves a lot of disk space, especially if used with the `--cleanup_work_dir` argument.

    ${c_yellow}Profile:${c_reset}
    You can merge different profiles for different setups, e.g.

        -profile local,docker
        -profile lsf,singularity
        -profile slurm,singularity

    -profile                 standard (local,docker) [default]

                             local
                             lsf
                             slurm

                             docker
                             singularity
                             conda
                             mamba

                             gcloud (use this as template for your own GCP setup)
                             ${c_reset}
    """.stripIndent()
}

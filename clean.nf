#!/usr/bin/env nextflow

if( !nextflow.version.matches('20.01+') ) {
    println "This workflow requires Nextflow version 20.01 or greater -- You are running version $nextflow.version"
    exit 1
}

nextflow.enable.dsl=2

/*
Nextflow -- Decontamination Pipeline
Author: marie.lataretu@uni-jena.de
Author: hoelzer.martin@gmail.com
*/

// Parameters sanity checking

Set valid_params = ['max_cores', 'cores', 'max_memory', 'memory', 'profile', 'help', 'nano', 'illumina', 'illumina_single_end', 'fasta', 'list', 'host', 'own', 'control', 'rm_rrna', 'bbduk', 'bbduk_kmer', 'reads_rna', 'output', 'nf_runinfo_dir', 'databases', 'condaCacheDir', 'singularityCacheDir', 'singularityCacheDir', 'cloudProcess', 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'] // don't ask me why there is also 'conda-cache-dir', 'singularity-cache-dir', 'cloud-process'
def parameter_diff = params.keySet() - valid_params
if (parameter_diff.size() != 0){
    exit 1, "ERROR: Parameter(s) $parameter_diff is/are not valid in the pipeline!\n"
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
Set hosts = ['hsa', 'mmu', 'cli', 'csa', 'gga', 'eco']

if (params.profile) { exit 1, "--profile is wrong, use -profile" }
if (params.nano == '' &&  params.illumina == '' && params.fasta == '' && params.illumina_single_end == '' ) { exit 1, "Read files missing, use [--nano] or [--illumina] or [--fasta]"}
if (params.control) { for( String ctr : params.control.split(',') ) if ( ! (ctr in controls ) ) { exit 1, "Wrong control defined (" + ctr + "), use one of these: " + controls } }
if (params.nano && params.control && 'dcs' in params.control.split(',') && 'eno' in params.control.split(',')) { exit 1, "Please choose either eno (for ONT dRNA-Seq) or dcs (for ONT DNA-Seq)." }
if (params.host) { for( String hst : params.host.split(',') ) if ( ! (hst in hosts ) ) { exit 1, "Wrong host defined (" + hst + "), use one of these: " + hosts } }
if (!params.host && !params.own && !params.control && !params.rm_rrna) { exit 1, "Please provide a control (--control), a host tag (--host), a FASTA file (--own) or set --rm_rrna for rRNA removal for the clean up."}

/************************** 
* INPUT CHANNELS 
**************************/

// nanopore reads input & --list support
if (params.nano && params.list) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] } 
} else if (params.nano) { nano_input_ch = Channel
    .fromPath( params.nano, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
}

// illumina paired-end reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
} else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
}

// illumina single-end reads input & --list support
if (params.illumina_single_end && params.list) { illumina_single_end_input_ch = Channel
  .fromPath( params.illumina_single_end, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
} else if (params.illumina_single_end) { illumina_single_end_input_ch = Channel
  .fromPath( params.illumina_single_end, checkIfExists: true )
  .map { file -> tuple(file.simpleName, file) } 
}

// assembly fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
} else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
}

// load control fasta sequence
if (params.control) {
  if ( 'phix' in params.control.split(',') ) {
    illuminaControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/phix.fa.gz' , checkIfExists: true )
  } else { illuminaControlFastaChannel = Channel.empty() }
  if ( 'dcs' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/dcs.fa.gz' , checkIfExists: true )
  } else if ( 'eno' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath( workflow.projectDir + '/data/controls/eno.fa.gz' , checkIfExists: true )
  } else { nanoControlFastaChannel = Channel.empty() }
} else {
  nanoControlFastaChannel = Channel.empty()
  illuminaControlFastaChannel = Channel.empty()
}

// load rRNA DB
if (params.rm_rrna){
  rRNAChannel = Channel.fromPath( workflow.projectDir + '/data/rRNA/*.fasta.gz', checkIfExists: true )
} else{
  rRNAChannel = Channel.empty()
}

if (params.host) {
  hostNameChannel = Channel.from( params.host ).splitCsv().flatten()
}

// user defined fasta sequence
if (params.own) {
  //TODO funzt auch mit * und so was?
  ownFastaChannel = Channel.fromPath( params.own, checkIfExists: true).splitCsv().flatten().map{ it -> file( it, checkIfExists: true ) }
}

/************************** 
* MODULES
**************************/

/* Comment section: */

include { download_host; check_own; concat_contamination } from './modules/get_host'

include { minimap2_fasta; minimap2_nano; minimap2_illumina } from './modules/minimap2'
include { bbduk } from './modules/bbmap'

include { rename_reads; restore_reads; get_number_of_reads; minimap2Stats; bbdukStats; writeLog } from './modules/utils'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow prepare_host {
  main:
    if (params.host) {
      if (params.cloudProcess) {
        host_preload = file("${params.databases}/hosts/${params.host}.fa.gz")
        if (host_preload.exists()) {
          host = Channel.fromPath(host_preload)
        } else {
          download_host(hostNameChannel)
          host = download_host.out
        }
      } else {
        download_host(hostNameChannel)
        host = download_host.out
      }
    }
    else {
      host = Channel.empty()
    }
    if (params.own) {
      check_own(ownFastaChannel)
      checkedOwn = check_own.out
    }
    else {
      checkedOwn = Channel.empty()
    }
  emit:
    host = host
    checkedOwn = checkedOwn
}

/************************** 
* SUB WORKFLOWS
**************************/

/* Comment section: */

workflow clean_fasta {
  take: 
    fasta_input_ch
    host
    checkedOwn
    rRNAChannel

  main:
    contamination = host.collect()
      .mix(illuminaControlFastaChannel)
      .mix(nanoControlFastaChannel)
      .mix(checkedOwn)
      .mix(rRNAChannel).collect()
    concat_contamination( contamination )
    minimap2_fasta(fasta_input_ch, concat_contamination.out)
    writeLog(fasta_input_ch.map{ it -> it[0] }, 'minimap2', fasta_input_ch.map{ it -> it[1] }, contamination)
    minimap2Stats(minimap2_fasta.out.name, minimap2_fasta.out.totalcontigs, minimap2_fasta.out.idxstats)
} 

workflow clean_nano {
  take: 
    nano_input_ch
    host
    checkedOwn
    rRNAChannel

  main:
    if (params.nano && params.illumina) {
      contamination = host.collect()
        .mix(nanoControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    } else {
      contamination = host.collect()
        .mix(nanoControlFastaChannel)
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    }
    rename_reads(nano_input_ch, 'single')
    minimap2_nano(rename_reads.out, concat_contamination.out)
    writeLog(nano_input_ch.map{ it -> it[0] }, 'minimap2', nano_input_ch.map{ it -> it[1] }, contamination)
    get_number_of_reads(rename_reads.out, 'single')
    minimap2Stats(minimap2_nano.out.idxstats.join(get_number_of_reads.out))
    restore_reads(minimap2_nano.out.cleaned_reads.concat(minimap2_nano.out.contaminated_reads), 'single', 'minimap2')
} 

workflow clean_illumina {
  take: 
    illumina_input_ch
    host
    checkedOwn
    rRNAChannel

  main:
    if (params.nano && params.illumina) {
      contamination = host.collect()
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    } else {
      contamination = host.collect()
        .mix(nanoControlFastaChannel)
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    }
    rename_reads(illumina_input_ch, 'paired')
    if (params.bbduk){
      bbduk(rename_reads.out, concat_contamination.out, 'paired')
      writeLog(illumina_input_ch.map{ it -> it[0] }, 'bbduk', illumina_input_ch.map{ it -> it[1] }, contamination)
      bbdukStats(bbduk.out.name, bbduk.out.stats)
      restore_reads(bbduk.out.cleaned_reads.concat(bbduk.out.contaminated_reads), 'paired', 'bbduk')
    } else {
      minimap2_illumina(rename_reads.out, concat_contamination.out, 'paired')
      writeLog(illumina_input_ch.map{ it -> it[0] }, 'minimap2', illumina_input_ch.map{ it -> it[1] }, contamination)
      get_number_of_reads(rename_reads.out, 'paired')
      minimap2Stats(minimap2_illumina.out.idxstats.join(get_number_of_reads.out))
      restore_reads(minimap2_illumina.out.cleaned_reads.concat(minimap2_illumina.out.contaminated_reads), 'paired', 'minimap2')
    }
} 

workflow clean_illumina_single {
  take:
    illumina_single_end_input_ch
    host
    checkedOwn
    rRNAChannel

  main:
    if (params.nano && params.illumina) {
      contamination = host.collect()
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    } else {
      contamination = host.collect()
        .mix(nanoControlFastaChannel)
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn)
        .mix(rRNAChannel).collect()
      concat_contamination( contamination )
    }
    rename_reads(illumina_single_end_input_ch, 'single')
    if (params.bbduk){
      bbduk(rename_reads.out, concat_contamination.out, 'single')
      writeLog(illumina_single_end_input_ch.map{ it -> it[0] }, 'bbduk', illumina_single_end_input_ch.map{ it -> it[1] }, contamination)
      bbdukStats(illumina_single_end_input_ch.map{ it -> it[0] }, bbduk.out.stats)
      restore_reads(bbduk.out.cleaned_reads.concat(bbduk.out.contaminated_reads), 'single', 'bbduk')
    } else {
      minimap2_illumina(rename_reads.out, concat_contamination.out, 'single')
      writeLog(illumina_single_end_input_ch.map{ it -> it[0] }, 'minimap2', illumina_single_end_input_ch.map{ it -> it[1] }, contamination)
      get_number_of_reads(rename_reads.out, 'single')
      minimap2Stats(minimap2_illumina.out.idxstats.join(get_number_of_reads.out))
      restore_reads(minimap2_illumina.out.cleaned_reads.concat(minimap2_illumina.out.contaminated_reads), 'single', 'minimap2')
    }
} 

/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
  prepare_host()

  if (params.fasta) {
    clean_fasta(fasta_input_ch, prepare_host.out.host, prepare_host.out.checkedOwn, rRNAChannel)
  }

  if (params.nano) { 
    clean_nano(nano_input_ch, prepare_host.out.host, prepare_host.out.checkedOwn, rRNAChannel)
  }

  if (params.illumina) { 
    clean_illumina(illumina_input_ch, prepare_host.out.host, prepare_host.out.checkedOwn, rRNAChannel)
  }

  if (params.illumina_single_end) { 
    clean_illumina_single(illumina_single_end_input_ch, prepare_host.out.host, prepare_host.out.checkedOwn, rRNAChannel)
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
    nextflow run clean.nf --nano '*/*.fastq' --host eco --control dcs 
    or
    nextflow run clean.nf --illumina '*/*.R{1,2}.fastq' --own some_host.fasta --bbduk 
    or
    nextflow run clean.nf --illumina 'test/illumina*.R{1,2}.fastq.gz' --nano data/nanopore.fastq.gz --fasta data/assembly.fasta --host eco --control phix

    ${c_yellow}Input:${c_reset}
    ${c_green}--nano ${c_reset}               '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green}--illumina ${c_reset}           '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green}--illumina_single_end${c_reset} '*.fastq.gz'                -> one sample per file
    ${c_green}--fasta ${c_reset}              '*.fasta.gz'                -> one sample per file
    ${c_dim} ...read above input from csv files:${c_reset} ${c_green}--list ${c_reset} 
                         ${c_dim}required format: name,path for --nano and --fasta; name,pathR1,pathR2 for --illumina; name,path for --illumina_single_end${c_reset}   

    ${c_yellow}Decontamination options:${c_reset}
    ${c_green}--host${c_reset}         comma separated list of reference genomes for decontamination, downloaded based on this parameter [default: $params.host]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - csa [NCBI: GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic]
                                        - gga [NCBI: Gallus_gallus.GRCg6a.dna.toplevel]
                                        - cli [NCBI: GCF_000337935.1_Cliv_1.0_genomic]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
    ${c_green}--control${c_reset}       comma separated list of common controls used in Illumina or Nanopore sequencing [default: $params.control]
                                        ${c_dim}Currently supported are:
                                        - phix [Illumina: enterobacteria_phage_phix174_sensu_lato_uid14015, NC_001422]
                                        - dcs [ONT DNA-Seq: a positive control (3.6 kb standard amplicon mapping the 3' end of the Lambda genome)]
                                        - eno [ONT RNA-Seq: a positive control (yeast ENO2 Enolase II of strain S288C, YHR174W)]${c_reset}
    ${c_green}--own ${c_reset}          use your own FASTA sequences (comma separated list of files) for decontamination, e.g. host.fasta.gz,spike.fasta [default: $params.own]
    ${c_green}--rm_rrna ${c_reset}      clean your data from rRNA [default: $params.rm_rrna]
    ${c_green}--bbduk${c_reset}         add this flag to use bbduk instead of minimap2 for decontamination of short reads [default: $params.bbduk]
    ${c_green}--bbduk_kmer${c_reset}    set kmer for bbduk [default: $params.bbduk_kmer]
    ${c_green}--reads_rna${c_reset}           add this flag for noisy direct RNA-Seq Nanopore data [default: $params.reads_rna]

    ${c_yellow}Compute options:${c_reset}
    --cores             max cores per process for local use [default $params.cores]
    --max_cores         max cores used on the machine for local use [default $params.max_cores]
    --memory            max memory for local use, enter in this format '8.GB' [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}Computing:${c_reset}
    In particular for execution of the workflow on a HPC (LSF, SLURM) adjust the following parameters:
    --databases             defines the path where databases are stored [default: $params.databases]
    --condaCacheDir         defines the path where environments (conda) are cached [default: $params.condaCacheDir]
    --singularityCacheDir   defines the path where images (singularity) are cached [default: $params.singularityCacheDir] 

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

                             ebi (lsf,singularity; preconfigured for the EBI cluster)
                             yoda (lsf,singularity; preconfigured for the EBI YODA cluster)
                             ara (slurm,conda; preconfigured for the ARA cluster)
                             gcloud (use this as template for your own GCP setup)
                             ${c_reset}
    """.stripIndent()
}

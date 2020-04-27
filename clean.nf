#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- Decontamination Pipeline
Author: marie.lataretu@uni-jena.de
Author: hoelzer.martin@gmail.com
*/

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
println "Workdir location:"
println "  $workflow.workDir\u001B[0m"
println " "
if (workflow.profile == 'standard') {
println "\033[2mCPUs to use: $params.cores"
println "Output dir name: $params.output\u001B[0m"
println " "}

if( !nextflow.version.matches('20.01+') ) {
    println "This workflow requires Nextflow version 20.01 or greater -- You are running version $nextflow.version"
    exit 1
}

Set controls = ['phix', 'dcs', 'eno']
Set hosts = ['hsa', 'mmu', 'cli', 'csa', 'gga', 'eco']

if (params.profile) { exit 1, "--profile is wrong, use -profile" }
if (params.nano == '' &&  params.illumina == '' && params.fasta == '' ) { exit 1, "Read files missing, use [--nano] or [--illumina] or [--fasta]"}
if (params.control) { for( String ctr : params.control.split(',') ) if ( ! (ctr in controls ) ) { exit 1, "Wrong control defined (" + ctr + "), use one of these: " + controls } }
if (params.nano && params.control && 'dcs' in params.control.split(',') && 'eno' in params.control.split(',')) { exit 1, "Please choose either eno (for ONT dRNA-Seq) or dcs (for ONT DNA-Seq)." }
if (params.host) { for( String hst : params.host.split(',') ) if ( ! (hst in hosts ) ) { exit 1, "Wrong host defined (" + hst + "), use one of these: " + hosts } }
if (!params.host && !params.own && !params.control) { exit 1, "Please provide a control (--control), a host tag (--host) or a FASTA file (--own) for the clean up."}

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

// illumina reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
} else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
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
    illuminaControlFastaChannel = Channel.fromPath(params.controldir + '/phix.fa.gz' , checkIfExists: true )
  } else { illuminaControlFastaChannel = Channel.empty() }
  if ( 'dcs' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath(params.controldir + '/dcs.fa.gz' , checkIfExists: true )
  } else if ( 'eno' in params.control.split(',') ) {
    nanoControlFastaChannel = Channel.fromPath(params.controldir + '/eno.fa.gz' , checkIfExists: true )
  } else { nanoControlFastaChannel = Channel.empty() }
}

if (params.host) {
  hostNameChannel = Channel.from( params.host ).splitCsv().flatten()
}

// user defined fasta sequence
if (params.own) {
  //TODO funzt auch mit * und so was?
  ownFastaChannel = Channel.from( params.own ).splitCsv().flatten().map{ it -> file( it, checkIfExists: true ) }
}

/************************** 
* MODULES
**************************/

/* Comment section: */

include {download_host; check_own; concat_contamination} from './modules/get_host'

include {minimap2_fasta; minimap2_nano; minimap2_illumina} from './modules/minimap2'
include {bbduk; bbdukStats} from './modules/bbmap'

/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow prepare {
  main:
    if (params.host) {
      download_host(hostNameChannel)
      host = download_host.out
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

// TODO only selective control when nano and illumina

workflow clean_fasta {
  take: 
    fasta_input_ch
    host
    checkedOwn

  main:
    concat_contamination(
      host.collect()
      .mix(illuminaControlFastaChannel)
      .mix(nanoControlFastaChannel)
      .mix(checkedOwn).collect())
    minimap2_fasta(fasta_input_ch, concat_contamination.out)
  
} 

workflow clean_nano {
  take: 
    nano_input_ch
    host
    checkedOwn

  main:
    if (params.nano && params.illumina) {
      concat_contamination(
        host.collect()
        .mix(nanoControlFastaChannel)
        .mix(checkedOwn).collect())
    } else {
      concat_contamination(
        host.collect()
        .mix(nanoControlFastaChannel)
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn).collect())
    }
    minimap2_nano(nano_input_ch, concat_contamination.out)
} 

workflow clean_illumina {
  take: 
    illumina_input_ch
    host
    checkedOwn

  main:
    if (params.nano && params.illumina) {
      concat_contamination(host.collect()
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn).collect())
    } else {
      concat_contamination(
        host.collect()
        .mix(nanoControlFastaChannel)
        .mix(illuminaControlFastaChannel)
        .mix(checkedOwn).collect())
    }
    if (params.bbduk){
      bbduk(illumina_input_ch, concat_contamination.out)
      bbdukStats(bbduk.out.name, bbduk.out.stats)
    } else {
      minimap2_illumina(illumina_input_ch, concat_contamination.out)
    }
} 


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
  prepare()

  if (params.fasta) {
    clean_fasta(fasta_input_ch, prepare.out.host, prepare.out.checkedOwn)
  }

  if (params.nano) { 
    clean_nano(nano_input_ch, prepare.out.host, prepare.out.checkedOwn)
  }

  if (params.illumina) { 
    clean_illumina(illumina_input_ch, prepare.out.host, prepare.out.checkedOwn)
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
    nextflow run clean.nf --illumina 'data/illumina*.R{1,2}.fastq.gz' --nano data/nanopore.fastq.gz --fasta data/assembly.fasta --host eco --control phix

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green} --fasta ${c_reset}           '*.fasta.gz'                -> one sample per file
    ${c_dim} ...read above input from csv files:${c_reset} ${c_green}--list ${c_reset} 
                         ${c_dim}required format: name,path for --nano and --fasta; name,pathR1,pathR2 for --illumina${c_reset}   

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
    ${c_green}--bbduk${c_reset}         add this flag to use bbduk instead of minimap2 for decontamination of short reads [default: $params.bbduk]
    ${c_green}--bbduk_kmer${c_reset}    set kmer for bbduk [default: $params.bbduk_kmer]
    ${c_green}--rna${c_reset}           add this flag for noisy direct RNA-Seq Nanopore data [default: $params.rna]

    ${c_yellow}Compute options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use, enter in this format '8.GB' [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 


    ${c_yellow}Profile:${c_reset}
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             ebi (HPC w/ LSF, singularity/docker, preconfigured for the EBI cluster)
                             gcloudMartin (googlegenomics and docker)
                             ${c_reset}
    """.stripIndent()
}

  
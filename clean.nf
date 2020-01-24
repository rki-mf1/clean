#!/usr/bin/env nextflow
nextflow.preview.dsl=2

/*
Nextflow -- Decontamination Pipeline
Author: hoelzer.martin@gmail.com
*/

/************************** 
* META & HELP MESSAGES 
**************************/

/* 
Comment section: First part is a terminal print for additional user information,
followed by some help statements (e.g. missing input) Second part is file
channel input. This allows via --list to alter the input of --nano & --illumina
to add csv instead. name,path   or name,pathR1,pathR2 in case of illumina 
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

if (params.profile) { exit 1, "--profile is WRONG use -profile" }
if (params.nano == '' &&  params.illumina == '' && params.fasta == '' ) { exit 1, "input missing, use [--nano] or [--illumina] or [--fasta]"}
if (params.control) {if (params.control != 'phix' &&  params.control != 'dcs' &&  params.control != 'eno') { exit 1, "wrong control defined, use [phix], [dcs] or [eno]"}}
if (!params.host && !params.own && !params.control) { exit 1, "please provide a control (--control), a host tag (--host) or a FASTA file (--own) for the clean up. A control can be combined with eiter --host or --own."}
if (params.host && params.own) {print "Attention, you provided a host via the --host flag (${params.host})\n and your own reference sequence via --own (${params.own}). Only your own file will be used.\n If you provided a control via --control (${params.control}) this will be added to your --own sequence.\n\n"}

/************************** 
* INPUT CHANNELS 
**************************/

// nanopore reads input & --list support
if (params.nano && params.list) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.nano) { nano_input_ch = Channel
    .fromPath( params.nano, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
}

// illumina reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  .view() }
  else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
}

// assembly fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
}

// user defined host genome fasta
host = false
if (params.own) {
  host = file(params.own, checkIfExists: true)
}

/************************** 
* MODULES
**************************/

/* Comment section: */

if (params.own) {
  include './modules/get_host' params(control: params.control, host: host.simpleName, own: params.own, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
} else {
  include './modules/get_host' params(control: params.control, host: params.host, own: params.own, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
}
if (params.own) {
  include './modules/build_bowtie2_index' params(control: params.control, host: host.simpleName, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
} else {
 include './modules/build_bowtie2_index' params(control: params.control, host: params.host, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
}

include './modules/minimap2' params(output: params.output)
include './modules/bowtie2' params(output: params.output, control: params.control)


/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow prepare_host {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { get_host(host); db = get_host.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (params.control) {
        if (params.host) {
          db_preload = file("${params.cloudDatabase}/hosts/${params.host}_${params.control}/${params.host}_${params.control}.fa.gz")
        } else {
          db_preload = file("${params.cloudDatabase}/hosts/${params.control}/${params.control}.fa.gz")
        }
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.host}/${params.host}.fa.gz")
      }
      if (db_preload.exists()) { db = db_preload }
      else  { get_host(host); db = get_host.out } 
    }
  emit: db
}

workflow bowtie2_index {
  get:
    genome

  main:
    // local storage via storeDir
    if (!params.cloudProcess) { build_bowtie2_index(genome); db = build_bowtie2_index.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (params.control) {
        db_preload = file("${params.cloudDatabase}/hosts/${params.host}_${params.control}/${params.host}_${params.control}/bt2")
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.host}/${params.host}/bt2")
      }
      if (db_preload.exists()) { db = db_preload }
      else  { build_bowtie2_index(genome); db = build_bowtie2_index.out } 
    }
  emit: db
}

/************************** 
* SUB WORKFLOWS
**************************/

/* Comment section: */

workflow clean_fasta {
  get: 
    fasta_input_ch
    db

  main:
    minimap2_fasta(fasta_input_ch, db)

} 

workflow clean_nano {
  get: 
    nano_input_ch
    db

  main:
    minimap2_nano(nano_input_ch, db)
} 

workflow clean_illumina {
  get: 
    illumina_input_ch
    db
    index

  main:
    if (params.bowtie){
      bowtie2_illumina(illumina_input_ch, db, index)
    } else {
      minimap2_illumina(illumina_input_ch, db)
    }
} 


/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
      prepare_host()
      host = prepare_host.out

      index = false
      if (params.bowtie) {
        bowtie2_index(host)
        index = bowtie2_index.out
      }

      if (params.fasta) { 
        clean_fasta(fasta_input_ch, host)
      }

      if (params.nano) { 
        clean_nano(nano_input_ch, host)
      }

      if (params.illumina) { 
        clean_illumina(illumina_input_ch, host, index)
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
    to a host but we recommend using the ${c_dim}--bowtie${c_reset} flag to switch to Bowtie2 to clean short-read data.  

    Use the ${c_dim}--host${c_reset} and ${c_dim}--control${c_reset} flag to download a host database or specify your ${c_dim}--own${c_reset} FASTA. 
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run clean.nf --nano '*/*.fastq' --host eco --control dcs 
    or
    nextflow run clean.nf --illumina '*/*.R{1,2}.fastq' --own some_host.fasta --bowtie 
    or
    nextflow run clean.nf --illumina 'data/illumina*.R{1,2}.fastq.gz' --nano data/nanopore.fastq.gz --fasta data/assembly.fasta --host eco --control phix

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green} --fasta ${c_reset}           '*.fasta.gz'                -> one sample per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Decontamination options:${c_reset}
    ${c_green}--host${c_reset}       reference genome for decontamination is downloaded based on this parameter [default: $params.host]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - csa [NCBI: GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic]
                                        - gga [NCBI: Gallus_gallus.GRCg6a.dna.toplevel]
                                        - cli [NCBI: GCF_000337935.1_Cliv_1.0_genomic]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
    ${c_green}--control${c_reset}       use one of these flags to remove commom controls used in Illumina or Nanopore sequencing [default: $params.control]
                                        ${c_dim}Currently supported are:
                                        - phix [Illumina: enterobacteria_phage_phix174_sensu_lato_uid14015, NC_001422]
                                        - dcs [ONT DNA-Seq: a positive control (3.6 kb standard amplicon mapping the 3' end of the Lambda genome)]${c_reset}
                                        - eno [ONT RNA-Seq: a positive control (yeast ENO2 Enolase II of strain S288C, YHR174W)]${c_reset}
    ${c_green}--own ${c_reset}          use your own FASTA sequence for decontamination, host.fasta.gz [default: $params.own]
    ${c_green}--bowtie${c_reset}        add this flag to use bowtie2 instead of minimap2 for decontamination of short reads [default: $params.bowtie]

    ${c_yellow}Compute options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.memory]
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

  
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

/************************** 
* INPUT CHANNELS 
**************************/

// nanopore reads input & --list support
if (params.nano && params.list) { nano_input_ch = Channel
  .fromPath( params.nano, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.nano) { nano_input_ch = Channel
    .fromPath( params.nano, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
}

// illumina reads input & --list support
if (params.illumina && params.list) { illumina_input_ch = Channel
  .fromPath( params.illumina, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  .view() }
  else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
  .view() 
}

// assembly fasta input & --list support
if (params.fasta && params.list) { fasta_input_ch = Channel
  .fromPath( params.fasta, checkIfExists: true )
  .splitCsv()
  .map { row -> ["${row[0]}", file("${row[1]}", checkIfExists: true)] }
  .view() }
  else if (params.fasta) { fasta_input_ch = Channel
    .fromPath( params.fasta, checkIfExists: true)
    .map { file -> tuple(file.simpleName, file) }
    .view()
}


/************************** 
* MODULES
**************************/

/* Comment section: */

include './modules/get_db' params(cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/module1' params(output: params.output, variable1: params.variable1)
include './modules/module2' params(output: params.output, variable2: params.variable2)


/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_db {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { example_db(); db = example_db.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      db_preload = file("${params.cloudDatabase}/test_db/Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa.gz")
      if (db_preload.exists()) { db = db_preload }
      else  { example_db(); db = example_db.out } 
    }
  emit: db
}


/************************** 
* SUB WORKFLOWS
**************************/

/* Comment section: */

workflow subworkflow_1 {
  get: 
    nano_input_ch
    db

  main:
    module2(module1(nano_input_ch, db))
} 

/************************** 
* WORKFLOW ENTRY POINT
**************************/

/* Comment section: */

workflow {
      download_db()
      db = download_db.out
      if (params.nano && !params.illumina) { 
        subworkflow_1(nano_input_ch, db)
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
    
    Workflow: Template
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run wf_template --nano '*/*.fastq' 

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green} --fasta ${c_reset}           '*.fasta.gz'                -> one sample per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_yellow}Parameters:${c_reset}
    --variable1             a variable [default: $params.variable1]
    --variable2             a variable [default: $params.variable2]

    ${c_dim}Nextflow options:
    -with-report rep.html    cpu / ram usage (may cause errors)
    -with-dag chart.html     generates a flowchart for the process tree
    -with-timeline time.html timeline (may cause errors)

    ${c_yellow}LSF computing:${c_reset}
    For execution of the workflow on a HPC with LSF adjust the following parameters:
    --databases         defines the path where databases are stored [default: $params.cloudDatabase]
    --workdir           defines the path where nextflow writes tmp files [default: $params.workdir]
    --cachedir          defines the path where images (singularity) are cached [default: $params.cachedir] 


    Profile:
    -profile                 standard (local, pure docker) [default]
                             conda (mixes conda and docker)
                             lsf (HPC w/ LSF, singularity/docker)
                             ebi (HPC w/ LSF, singularity/docker, preconfigured for the EBI cluster)
                             gcloudMartin (googlegenomics and docker)
                             ${c_reset}
    """.stripIndent()
}

  
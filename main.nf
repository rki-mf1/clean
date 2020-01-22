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
  .map { row -> [row[0], file("${row[1]}", checkIfExists: true)] }
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
  .map { row -> [row[0], [file("${row[1]}", checkIfExists: true), file("${row[2]}", checkIfExists: true)]] }
  .view() }
  else if (params.illumina) { illumina_input_ch = Channel
  .fromFilePairs( params.illumina , checkIfExists: true )
  .view() 
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
    .view()
}


/************************** 
* MODULES
**************************/

/* Comment section: */

include './modules/get_host' params(phix: params.phix, species: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)
include './modules/build_bowtie2_index' params(phix: params.phix, species: params.species, cloudProcess: params.cloudProcess, cloudDatabase: params.cloudDatabase)

include './modules/minimap2' params(output: params.output)
include './modules/bowtie2' params(output: params.output, phix: params.phix)


/************************** 
* DATABASES
**************************/

/* Comment section:
The Database Section is designed to "auto-get" pre prepared databases.
It is written for local use and cloud use via params.cloudProcess.
*/

workflow download_genomes {
  main:
    // local storage via storeDir
    if (!params.cloudProcess) { get_host(); db = get_host.out }
    // cloud storage via db_preload.exists()
    if (params.cloudProcess) {
      if (params.phix) {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}_phix/${params.species}_phix.fa.gz")
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}/${params.species}.fa.gz")
      }
      if (db_preload.exists()) { db = db_preload }
      else  { get_host(); db = get_host.out } 
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
      if (params.phix) {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}_phix/${params.species}_phix/bt2")
      } else {
        db_preload = file("${params.cloudDatabase}/hosts/${params.species}/${params.species}/bt2")
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
      download_genomes()
      genomes = download_genomes.out

      index = false
      if (params.bowtie) {
        bowtie2_index(genomes)
        index = bowtie2_index.out
      }

      if (params.fasta && !params.nano && !params.illumina) { 
        clean_fasta(fasta_input_ch, genomes)
      }

      if (!params.fasta && params.nano && !params.illumina) { 
        clean_nano(nano_input_ch, genomes)
      }

      if (!params.fasta && !params.nano && params.illumina) { 
        clean_illumina(illumina_input_ch, genomes, index)
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
    
    ${c_yellow}Usage example:${c_reset}
    nextflow run main.nf --nano '*/*.fastq' --species mmu --phix 
    or
    nextflow run main.nf --illumina '*/*.R{1,2}.fastq' --species eco --bowtie 

    ${c_yellow}Input:${c_reset}
    ${c_green} --nano ${c_reset}            '*.fasta' or '*.fastq.gz'   -> one sample per file
    ${c_green} --illumina ${c_reset}        '*.R{1,2}.fastq.gz'         -> file pairs
    ${c_green} --fasta ${c_reset}           '*.fasta.gz'                -> one sample per file
    ${c_dim}  ..change above input to csv:${c_reset} ${c_green}--list ${c_reset}            

    ${c_yellow}Options:${c_reset}
    --cores             max cores for local use [default: $params.cores]
    --memory            max memory for local use [default: $params.memory]
    --output            name of the result folder [default: $params.output]

    ${c_green}--species${c_reset}       reference genome for decontamination is selected based on this parameter [default: $params.species]
                                        ${c_dim}Currently supported are:
                                        - hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
                                        - mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
                                        - cli [NCBI: GCF_000337935.1_Cliv_1.0_genomic]
                                        - eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
    ${c_green}--phix${c_reset}          add this flag to download and add phiX genome for decontamination [default: $params.phix]
    ${c_green}--bowtie${c_reset}        add this flag to use bowtie2 instead of minimap2 for decontamination of short reads [default: $params.bowtie]


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

  
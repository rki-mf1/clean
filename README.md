# Clean your data!
A decontamination workflow for short reads, long reads and assemblies. 

![](https://img.shields.io/badge/nextflow-19.10.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)

Maintainer: Martin

Email: hoelzer.martin@gmail.com

# Objective
Sequencing data is often contaminated with DNA or RNA from other species. These, normally unwanted, material occurs for biological 
reasons or can be also spiked in as a control. For example, this is often the case for Illumina data ([phiX phage](https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/1944-3277-10-18)) or Oxford Nanopore
Technologies ([DNA CS (DCS)](https://assets.ctfassets.net/hkzaxo8a05x5/2IX56YmF5ug0kAQYoAg2Uk/159523e326b1b791e3b842c4791420a6/DNA_CS.txt)) 

# What this workflow does for you
With this workflow you can clean your Illumina, Nanopore or any FASTA-formated sequence date. The output are the clean and as contaminated identified sequences. 
Per default [minimap2](https://github.com/lh3/minimap2) is used for aligning your sequences to a host but I recommend using the [Bowtie2](https://github.com/BenLangmead/bowtie2) to clean short-read data (_--bowtie_).  

You can simply specify provided hosts and controls for the cleanup or use your own FASTA.    

# How does it work?

## Installation

* runs with the workflow manager `nextflow` using `docker` or `conda`
* this means all programs are automatically pulled via docker or conda
* only `docker` and `nextflow` need to be installed

### Easy 
If you dont have experience with bioinformatic tools just copy the commands into your terminal to set everything up:
```bash
sudo apt-get update
sudo apt install -y default-jre
curl -s https://get.nextflow.io | bash 
sudo mv nextflow /bin/
sudo apt-get install -y docker-ce docker-ce-cli containerd.io
sudo usermod -a -G docker $USER
```

* restart your computer
* try out the installation by entering the following

```bash
nextflow run hoelzer/clean --nano ~/.nextflow/assets/hoelzer/clean/data/nanopore.fastq.gz --host eco
```

### Experienced

**Dependencies**

>   * docker (add docker to your Usergroup, so no sudo is needed)
>   * nextflow + java runtime 
>   * git (should be already installed)
>   * wget (should be already installed)
>   * tar (should be already installed)

* Docker installation [here](https://docs.docker.com/v17.09/engine/installation/linux/docker-ce/ubuntu/#install-docker-ce)
* Nextflow installation [here](https://www.nextflow.io/)
* move or add the nextflow executable to a bin path
* add docker to your User group via `sudo usermod -a -G docker $USER`

# Execution examples

Get help:
```bash
nextflow run hoelzer/clean --help
```

Clean Nanopore data by filtering against a combined reference of the _E. coli_ genome and the Nanopore DNA CS spike-in.  
```bash
nextflow run hoelzer/clean --nano '*/*.fastq.gz' --host eco --control dcs 
```

Clean Illumina paired-end data against your own reference FASTA using Bowtie2 instead of minimap2. 
```bash
nextflow run hoelzer/clean --illumina '*/*.R{1,2}.fastq' --own some_host.fasta --bowtie 
```

Clean some Illumina, Nanopore, and assembly files against the mouse and phiX genomes.  
```bash
nextflow run hoelzer/clean --illumina 'data/illumina*.R{1,2}.fastq.gz' --nano data/nanopore.fastq.gz --fasta data/assembly.fasta --host mmu --control phix
```

# Supported species
Currently supported are:
* hsa | _Homo sapiens_ | [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
* mmu | _Mus musculus_ | [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
* csa | _Chlorocebus sabeus_ | [NCBI: GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic]
* gga | _Gallus gallus_ | [NCBI: Gallus_gallus.GRCg6a.dna.toplevel]
* cli | _Columba livia_ | [NCBI: GCF_000337935.1_Cliv_1.0_genomic]
* eco | _Escherichia coli_ | [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
... for reasons. More can be easily added! Just write me, add an issue or make a pull request. 

# Flowchart
![chart](figures/dag.png)
# Clean your data!
A decontamination workflow for short reads, long reads and assemblies. 

![](https://img.shields.io/badge/nextflow-19.10.0-brightgreen)
![](https://img.shields.io/badge/uses-docker-blue.svg)

Maintainer: Martin

Email: hoelzer.martin@gmail.com


# Input examples

* **one** .fastq file per sample: `--nano 'sample1.fastq'`
* paired end illumina: `--illumina 'S_41_17_Cf*.R{1,2}.fastq.gz'`
* **one** .fasta file per sample: `--fasta 'sample1.fasta'`

# Execution example

````
nextflow run main.nf
````

# Supported species
Currently supported are:
* hsa [Ensembl: Homo_sapiens.GRCh38.dna.primary_assembly]
* mmu [Ensembl: Mus_musculus.GRCm38.dna.primary_assembly]
* csa [NCBI: GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic]
* gga [NCBI: Gallus_gallus.GRCg6a.dna.toplevel]
* cli [NCBI: GCF_000337935.1_Cliv_1.0_genomic]
* eco [Ensembl: Escherichia_coli_k_12.ASM80076v1.dna.toplevel]${c_reset}
... for reasons. More can be easily added!

# Flowchart
![chart](figures/dag.png)
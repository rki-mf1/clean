/*Comment section: */

process bowtie2_illumina {
  label 'bowtie2'
  publishDir "${params.output}/${name}/bowtie", mode: 'copy', pattern: "*.gz" 
  publishDir "${params.output}/${name}/bowtie", mode: 'copy', pattern: "log.txt" 

  input: 
    tuple val(name), file(reads)
    file(genome)
    file(index)

  output:
    tuple file("*.gz"), file('log.txt')

  script:
    """
    #cp bt2/*.bt2 .
    bowtie2 -x ${genome.simpleName} -1 ${reads[0]} -2 ${reads[1]} -S ${name}.sam

    samtools fastq -F 2 -1 ${name}.clean.R1.fastq -2 ${name}.clean.R2.fastq ${name}.sam
    samtools fastq -f 2 -1 ${name}.contamination.R1.fastq -2 ${name}.contamination.R2.fastq ${name}.sam

    rm *.gz *.bt2
    gzip ${name}.clean.R1.fastq
    gzip ${name}.clean.R2.fastq
    gzip ${name}.contamination.R1.fastq
    gzip ${name}.contamination.R2.fastq

    touch log.txt
    cat <<EOF >> log.txt
Input:\t${reads[0]}, ${reads[1]} 
Host:\t${genome}

Clean:\t\t${params.output}/${name}/minimap2/${name}.clean.R1.fastq.gz
\t\t${params.output}/${name}/minimap2/${name}.clean.R2.fastq.gz
Contaminated:\t${params.output}/${name}/minimap2/${name}.contamination.R1.fastq.gz
\t\t${params.output}/${name}/minimap2/${name}.contamination.R2.fastq.gz

# Stay clean!
EOF
    """
}

process bowtie2_illumina_f12 {
  label 'bowtie2'
  publishDir "${params.output}/${name}/bowtie_f12", mode: 'copy', pattern: "*.gz" 
  publishDir "${params.output}/${name}/bowtie", mode: 'copy', pattern: "log.txt" 

  input: 
    tuple val(name), file(reads)
    file(genome)
    file(index)

  output:
    tuple file("*.gz"), file('log.txt')

  script:
    """
    #cp bt2/*.bt2 .
    bowtie2 -x ${genome.simpleName} -1 ${reads[0]} -2 ${reads[1]} -S ${name}.sam

    samtools fastq -f 12 -F 256 -1 ${name}.clean.R1.fastq -2 ${name}.clean.R2.fastq ${name}.sam

    rm *.gz *.bt2
    gzip ${name}.clean.R1.fastq
    gzip ${name}.clean.R2.fastq

    touch log.txt
    cat <<EOF >> log.txt
Input:\t${reads[0]}, ${reads[1]} 
Host:\t${genome}

Clean:\t\t${params.output}/${name}/minimap2/${name}.clean.R1.fastq.gz
\t\t${params.output}/${name}/minimap2/${name}.clean.R2.fastq.gz
Contaminated:\t${params.output}/${name}/minimap2/${name}.contamination.R1.fastq.gz
\t\t${params.output}/${name}/minimap2/${name}.contamination.R2.fastq.gz

# Stay clean!
EOF
    """
}

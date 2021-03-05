/*Comment section: */

process minimap2_fasta {
  label 'minimap2'

  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.gz"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam"
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam.bai"

  input: 
    tuple val(name), path(fasta)
    path db

  output:
    val name, emit: name
    path '*.gz'
    path '*.contamination.sorted.bam'
    path '*.contamination.sorted.bam.bai'
    path 'idxstats.tsv', emit: idxstats
    env TOTALCONTIGS, emit: totalcontigs

  script:
  """
  if [[ ${fasta} =~ \\.gz\$ ]]; then
    TOTALCONTIGS=\$(zgrep '^>' ${fasta} | wc -l)
  else
    TOTALCONTIGS=\$(grep '^>' ${fasta} | wc -l)
  fi

  minimap2 -ax asm5 -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${fasta}

  samtools fasta -f 4 -0 ${name}.clean.fasta ${name}.sam
  samtools fasta -F 4 -0 ${name}.contamination.fasta ${name}.sam
  pigz -p ${task.cpus} ${name}.clean.fasta
  pigz -p ${task.cpus} ${name}.contamination.fasta

  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv

  rm -f ${name}.sam
  """
}

process minimap2_nano {
  label 'minimap2'
  
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam*"

  input: 
    tuple val(name), path(reads)
    path db

  output:
    tuple val(name), path ('idxstats.tsv'), emit: idxstats
    tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
    tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads
    path '*.contamination.sorted.bam*'

  script:
  """
  PARAMS="-ax map-ont"
  if [[ ${params.reads_rna} != 'false' ]]; then
    PARAMS="-ax splice -k14"
  fi

  minimap2 \$PARAMS -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads}

  samtools fastq -f 4 -0 ${reads.baseName}.clean.fastq ${name}.sam
  samtools fastq -F 4 -0 ${reads.baseName}.contamination.fastq ${name}.sam

  samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
  samtools index ${name}.contamination.sorted.bam
  samtools idxstats  ${name}.contamination.sorted.bam > idxstats.tsv
  """
}

process minimap2_illumina {
  label 'minimap2'

  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "*.contamination.sorted.bam*"

  input: 
    tuple val(name), path(reads)
    path db
    val mode

  output:
    tuple val(name), path ('idxstats.tsv'), emit: idxstats
    tuple val(name), val('clean'), path('*clean.fastq'), emit: cleaned_reads
    tuple val(name), val('contamination'), path('*contamination.fastq'), emit: contaminated_reads
    path '*.contamination.sorted.bam*'

  script:
  if ( mode == 'paired' ) {
    """
    minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads[0]} ${reads[1]}
    
    # Use samtools -F 2 to discard only reads mapped in proper pair:
    samtools fastq -F 2 -1 ${reads[0].baseName}.clean.fastq -2 ${reads[1].baseName}.clean.fastq ${name}.sam
    samtools fastq -f 2 -1 ${reads[0].baseName}.contamination.fastq -2 ${reads[1].baseName}.contamination.fastq ${name}.sam

    samtools view -b -f 2 -F 2048 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv
    """
  } else {
    """
    minimap2 -ax sr -N 5 --secondary=no -t ${task.cpus} -o ${name}.sam ${db} ${reads}
    
    samtools fastq -f 4 -0 ${reads.baseName}.clean.fastq ${name}.sam
    samtools fastq -F 4 -0 ${reads.baseName}.contamination.fastq ${name}.sam

    samtools view -b -F 2052 ${name}.sam | samtools sort -o ${name}.contamination.sorted.bam --threads ${task.cpus}
    samtools index ${name}.contamination.sorted.bam
    samtools idxstats ${name}.contamination.sorted.bam > idxstats.tsv
    """
  }
}

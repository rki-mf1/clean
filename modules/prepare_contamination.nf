process download_host {
  label 'minimap2'

  if (params.cloudProcess) {
    publishDir "${params.databases}/hosts", mode: params.publish_dir_mode, pattern: "*.fa.gz"
  }
  else {
    storeDir "${params.databases}/hosts"
  }

  input:
  val host

  output:
  path "${host}.fa.gz"

  script:
  """
  case $host in
    hsa)
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -O host-temp.fa.gz
      ;;
    mmu)
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz -O host-temp.fa.gz
      ;;
    cli)
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz -O host-temp.fa.gz
      ;;
    csa)
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/409/795/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz -O host-temp.fa.gz
      ;;
    gga)
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz -O host-temp.fa.gz
      ;;
    eco)
      wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz -O host-temp.fa.gz
      ;;
    *)
      echo "Unknown host ($host)."
      ;;
  esac

  zcat < host-temp.fa.gz | bgzip -@ ${task.cpus} -c > ${host}.fa.gz
  """
  stub:
  """
  touch ${host}.fa.gz
  """
}

process check_own {
  label 'seqkit'

  input:
  path fasta

  output:
  path '*.gz', includeInputs: true

  script:
  """
  seqkit seq ${fasta} -o ${fasta}.gz
  """
  stub:
  """
  touch ${fasta}.gz
  """
}

process concat_contamination {
  label 'seqkit'

  publishDir (
    path: "${params.output}/intermediate",
    mode: params.publish_dir_mode,
    pattern: "db.fa.gz",
    enabled: !params.no_intermediate,
    saveAs: { "host.fa.gz" }
  )
  publishDir (
    path: "${params.output}/intermediate",
    mode: params.publish_dir_mode,
    pattern: "db.fa.fai",
    enabled: !params.no_intermediate,
    saveAs: { "host.fa.fai" }
  )

  input:
  path fastas

  output:
  path 'db.fa.gz', emit: fa
  path 'db.fa.fai', emit: fai

  script:
  """
  # Combine input files, rename duplicate sequences (by id) if found, and compress
  seqkit seq ${fastas} | seqkit rename | bgzip -@ ${task.cpus} -c > db.fa.gz
  samtools faidx db.fa.gz --fai-idx db.fa.fai
  """
  stub:
  """
  touch db.fa.gz db.fa.fai
  """
}

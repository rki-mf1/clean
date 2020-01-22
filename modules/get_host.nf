/*Comment section: */

process get_host {
  label 'basics'
  if (params.cloudProcess) { 
    if (params.phix) {
      publishDir "${params.cloudDatabase}/hosts/${params.species}_phix", mode: 'copy', pattern: "*.fa.gz" 
    } else {
      publishDir "${params.cloudDatabase}/hosts/${params.species}", mode: 'copy', pattern: "*.fa.gz" 
    }
  }
  else { 
    if (params.phix) {
      storeDir "nextflow-autodownload-databases/hosts/${params.species}_phix" 
    } else {
      storeDir "nextflow-autodownload-databases/hosts/${params.species}" 
    }
  }  

  output:
      file("*.fa.gz")

  script:
    """
    if [ ${params.species} == 'hsa' ]; then
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'mmu' ]; then
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'cli' ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz
      mv *.fna.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'csa' ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/409/795/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz
      mv *.fna.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'gga' ]; then
      wget ftp://ftp.ensembl.org/pub/release-99/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ ${params.species} == 'eco' ]; then
      wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
      mv *.fa.gz ${params.species}.fa.gz
    fi
    if [ "${params.phix}" != "false" ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna
      gzip -f NC_001422.fna
      cat *.gz > ${params.species}_phix.fa.gz.1
      rm *.gz
      mv ${params.species}_phix.fa.gz.1 ${params.species}_phix.fa.gz
    fi
    """
}


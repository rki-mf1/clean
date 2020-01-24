
process get_host {
  label 'basics'
  if (params.cloudProcess) { 
    if (params.control) {
      if (params.host) {
        publishDir "${params.cloudDatabase}/hosts/${params.host}_${params.control}", mode: 'copy', pattern: "*.fa.gz" 
      } else {
        publishDir "${params.cloudDatabase}/hosts/${params.control}", mode: 'copy', pattern: "*.fa.gz" 
      }
    } else {
      publishDir "${params.cloudDatabase}/hosts/${params.host}", mode: 'copy', pattern: "*.fa.gz" 
    }
  }
  else { 
    if (params.control) {
      if (params.host) {
        storeDir "nextflow-autodownload-databases/hosts/${params.host}_${params.control}" 
      } else {
        storeDir "nextflow-autodownload-databases/hosts/${params.control}" 
      }
    } else {
      storeDir "nextflow-autodownload-databases/hosts/${params.host}" 
    }
  }  

  input:
    file(own_fasta)

  output:
      file("*.fa.gz")

  script:
    """

    if [ ${params.own} == 'false' ]; then
      ## DOWNLOAD GENOME
     if [ ${params.host} == 'hsa' ]; then
       wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
       mv *.fa.gz ${params.host}.fa.gz
     fi
     if [ ${params.host} == 'mmu' ]; then
       wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
       mv *.fa.gz ${params.host}.fa.gz
     fi
     if [ ${params.host} == 'cli' ]; then
       wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz
       mv *.fna.gz ${params.host}.fa.gz
     fi
      if [ ${params.host} == 'csa' ]; then
       wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/409/795/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz
       mv *.fna.gz ${params.host}.fa.gz
     fi
     if [ ${params.host} == 'gga' ]; then
       wget ftp://ftp.ensembl.org/pub/release-99/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
       mv *.fa.gz ${params.host}.fa.gz
     fi
     if [ ${params.host} == 'eco' ]; then
       wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
       mv *.fa.gz ${params.host}.fa.gz
     fi
    else
      ## USE USER GENOME
      if [[ ${params.own} =~ \\.gz\$ ]]; then
        cp \$(basename ${params.own}) ${params.host}.fa.gz
      else
        cp \$(basename ${params.own}) ${params.host}.fa
        gzip -f ${params.host}.fa
      fi
    fi



    if [ "${params.control}" == "phix" ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna
      gzip -f NC_001422.fna
      if [ ${params.host} == 'false' ]; then
        mv NC_001422.fna.gz ${params.control}.fa.gz
      else
        cat *.gz > ${params.host}_${params.control}.fa.gz.1
        rm *.gz
        mv ${params.host}_${params.control}.fa.gz.1 ${params.host}_${params.control}.fa.gz
      fi
    fi
    
    if [ "${params.control}" == "dcs" ]; then
      wget https://assets.ctfassets.net/hkzaxo8a05x5/2IX56YmF5ug0kAQYoAg2Uk/159523e326b1b791e3b842c4791420a6/DNA_CS.txt
      echo "" >> DNA_CS.txt
      mv DNA_CS.txt DNA_CS.fa
      gzip -f DNA_CS.fa
      if [ ${params.host} == 'false' ]; then
        mv DNA_CS.fa.gz ${params.control}.fa.gz
      else
        cat *.gz > ${params.host}_${params.control}.fa.gz.1
        rm *.gz
        mv ${params.host}_${params.control}.fa.gz.1 ${params.host}_${params.control}.fa.gz
      fi
    fi
    """
}

/*

*/

/*
    HOST=${params.host}
    if [ ${params.own} != 'false' ]; then
      HOST=\$(basename ${params.own})
      mv \$HOST ${} 
      gzip -f \$HOST
      ##TODO what if the file is gzipped!
    else
      if [ \$HOST == 'hsa' ]; then
        wget ftp://ftp.ensembl.org/pub/release-99/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        mv *.fa.gz \$HOST.fa.gz
     fi
      if [ \$HOST == 'mmu' ]; then
        wget ftp://ftp.ensembl.org/pub/release-99/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
        mv *.fa.gz \$HOST.fa.gz
      fi
      if [ \$HOST == 'cli' ]; then
       wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/337/935/GCF_000337935.1_Cliv_1.0/GCF_000337935.1_Cliv_1.0_genomic.fna.gz
       mv *.fna.gz \$HOST.fa.gz
     fi
      if [ \$HOST == 'csa' ]; then
        wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/409/795/GCF_000409795.2_Chlorocebus_sabeus_1.1/GCF_000409795.2_Chlorocebus_sabeus_1.1_genomic.fna.gz
        mv *.fna.gz \$HOST.fa.gz
      fi
     if [ \$HOST == 'gga' ]; then
        wget ftp://ftp.ensembl.org/pub/release-99/fasta/gallus_gallus/dna/Gallus_gallus.GRCg6a.dna.toplevel.fa.gz
        mv *.fa.gz \$HOST.fa.gz
      fi
      if [ \$HOST == 'eco' ]; then
        wget ftp://ftp.ensemblgenomes.org/pub/release-45/bacteria//fasta/bacteria_90_collection/escherichia_coli_k_12/dna/Escherichia_coli_k_12.ASM80076v1.dna.toplevel.fa.gz
       mv *.fa.gz \$HOST.fa.gz
      fi
    fi

    if [ "${params.control}" == "phix" ]; then
      wget ftp://ftp.ncbi.nlm.nih.gov/genomes/Viruses/enterobacteria_phage_phix174_sensu_lato_uid14015/NC_001422.fna
      gzip -f NC_001422.fna
      if [ \$HOST == 'false' ]; then
        mv NC_001422.fna.gz ${params.control}.fa.gz
      else
        cat *.gz > \$HOST_${params.control}.fa.gz.1
        rm *.gz
        mv \$HOST_${params.control}.fa.gz.1 \$HOST_${params.control}.fa.gz
      fi
    fi
    
    if [ "${params.control}" == "dcs" ]; then
      wget https://assets.ctfassets.net/hkzaxo8a05x5/2IX56YmF5ug0kAQYoAg2Uk/159523e326b1b791e3b842c4791420a6/DNA_CS.txt
      mv DNA_CS.txt DNA_CS.fa
      gzip -f DNA_CS.fa
      if [ \$HOST == 'false' ]; then
        mv DNA_CS.fa.gz ${params.control}.fa.gz
      else
        cat *.gz > \$HOST_${params.control}.fa.gz.1
        rm *.gz
        mv \$HOST_${params.control}.fa.gz.1 \$HOST_${params.control}.fa.gz
      fi
    fi

*/
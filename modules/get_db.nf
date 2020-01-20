/*Comment section: */

process example_db {
  label 'label'
  if (params.cloudProcess) { 
    publishDir "${params.cloudDatabase}/test_db/", mode: 'copy', pattern: "Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa" 
  }
  else { 
    storeDir "nextflow-autodownload-databases/test_db/" 
  }  

  output:
    file("Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa")

  script:
    """
    wget ftp://ftp.ensemblgenomes.org/pub/bacteria/release-45/fasta/bacteria_44_collection/chlamydia_gallinacea_08_1274_3/dna/Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa.gz
    gunzip Chlamydia_gallinacea_08_1274_3.ASM47102v2.dna.toplevel.fa.gz
    """
}


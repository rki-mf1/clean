/*Comment section: */

process build_bowtie2_index {
  label 'bowtie2'
  if (params.cloudProcess) { 
    if (params.phix) {
      publishDir "${params.cloudDatabase}/hosts/${params.species}_phix/bt2", mode: 'copy', pattern: "*.bt2" 
    } else {
      publishDir "${params.cloudDatabase}/hosts/${params.species}/bt2", mode: 'copy', pattern: "*.bt2" 
    }
  }
  else { 
    if (params.phix) {
      storeDir "nextflow-autodownload-databases/hosts/${params.species}_phix/bt2" 
    } else {
      storeDir "nextflow-autodownload-databases/hosts/${params.species}/bt2" 
    }
  }  

  input:
    file(genome)

  output:
    file("*.bt2")

  script:
    """
    bowtie2-build ${genome} ${genome.simpleName}
    #mkdir bt2
    #mv *.bt2 bt2/
    """
}


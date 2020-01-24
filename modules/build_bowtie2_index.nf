/*Comment section: */

process build_bowtie2_index {
  label 'bowtie2'
  if (params.cloudProcess) { 
    if (params.control) {
      if (params.host) {
        publishDir "${params.cloudDatabase}/hosts/${params.host}_${params.control}/bt2", mode: 'copy', pattern: "*.bt2"
      } else {
        publishDir "${params.cloudDatabase}/hosts/${params.control}/bt2", mode: 'copy', pattern: "*.bt2"
      } 
    } else {
      publishDir "${params.cloudDatabase}/hosts/${params.host}/bt2", mode: 'copy', pattern: "*.bt2" 
    }
  }
  else { 
    if (params.control) {
      if (params.host) {
        storeDir "nextflow-autodownload-databases/hosts/${params.host}_${params.control}/bt2" 
      } else {
        storeDir "nextflow-autodownload-databases/hosts/${params.control}/bt2" 
      }
    } else {
      storeDir "nextflow-autodownload-databases/hosts/${params.host}/bt2" 
    }
  }  

  input:
    file(genome)

  output:
    file("*.bt2")

  script:
    """
    bowtie2-build ${genome} ${genome.simpleName}
    """
}


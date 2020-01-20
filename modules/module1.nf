/*Comment section: */

process module1 {
  label 'label'  
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}.${params.variable1}.fasta"

  input:
    tuple val(name), file(read)
    file(db)

  output:
    tuple val(name), file("${name}.${params.variable1}.fasta")

  script:
    """
    cat ${read} ${db} > ${name}.${params.variable1}.fasta
    """
}



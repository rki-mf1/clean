/*Comment section: */

process module2 {
  label 'label'  
  publishDir "${params.output}/${name}/", mode: 'copy', pattern: "${name}.${params.variable2}.gz"
  
  input:
    tuple val(name), file(combined)
  
  output:
    tuple val(name), file("${name}.${params.variable2}.gz")
  
  script:
    """
    gzip -f ${combined}
    mv ${combined}.gz ${name}.${params.variable2}.gz
    """
}

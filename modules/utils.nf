process minimap2Stats {
  label 'minimap2'
  publishDir "${params.output}/${name}/minimap2", mode: 'copy', pattern: "stats.txt" 

  input:
  val name
  val totalreads
  path idxstats

  output:
  path "stats.txt"

  script:
  """
  MAPPEDSUM=\$(awk -F '\\t' '{sum += \$3} END {print sum}' idxstats.tsv)
  UNPROMAPPEDSUM=\$(awk -F '\\t' '/^[^*]/ {sum += \$4} END {print sum}' idxstats.tsv)

  MAP=\$(awk -v map=\$MAPPEDSUM -v unpromap=\$UNPROMAPPEDSUM -v tot=${totalreads} 'BEGIN {perc=(map-unpromap)/tot*100; print map-unpromap " ("  perc " %) reads were properly mapped; of these:"}')

  FA=\$(awk -v tot=${totalreads} -F '\\t' '/^[^*]/ {propmap=\$3-\$4; print "\\t\\t" propmap " (" propmap/tot*100  "%) reads aligned to " \$1}' idxstats.tsv)

  touch stats.txt
  cat <<EOF >> stats.txt
  ${totalreads} reads in total; of these:
  \t\$MAP
  \$FA
  EOF
  """
}

process bbdukStats {
  publishDir "${params.output}/${name}/bbduk", mode: 'copy', pattern: "stats.txt"

  input:
  val name
  path bbdukStats

  output:
  path "stats.txt"

  script:
  """
  TOTAL=\$(grep '#Total' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MNUM=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$2}')
  MPER=\$(grep '#Matched' ${bbdukStats} | awk -F '\\t' '{print \$3}')

  FA=\$(awk -F '\\t' '/^[^#]/ {print "\\t\\t"\$2" ("\$3") aligned to "\$1}' ${bbdukStats})

  touch stats.txt
  cat <<EOF >> stats.txt
  \$TOTAL reads in total; of these:
  \t\$MNUM (\$MPER) reads were properly mapped; of these:
  \$FA
  EOF
  """
}

process writeLog {
  publishDir "${params.output}/${name}/${tool}", mode: 'copy', pattern: "log.txt"

  input:
    val name
    val tool
    val reads
    val db

  output:
    path 'log.txt'
  
  script:
  """
  touch log.txt
  cat <<EOF >> log.txt
  Input reads:\t${reads}
  Contamination:\t${db}
  
  Statistics summary:\t${params.output}/${name}/${tool}/stats.txt
  EOF
  """

}
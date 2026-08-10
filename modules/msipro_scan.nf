process MSIPRO_SCAN {
  tag { ref_fa.baseName }
  container 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'

  // global reference outputs
  publishDir { "${params.outdir}/reference/msisensor" }, mode: 'copy'

  input:
  path ref_fa

  output:
  path "microsat_sites.tsv", emit: sites
  path ref_fa,               emit: fasta

  script:
  """
  set -euo pipefail

  msisensor-pro scan \
    -d ${ref_fa} \
    -o microsat_sites.tsv \
    -p ${task.cpus}
  """
}

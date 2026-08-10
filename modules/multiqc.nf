process MULTIQC {
  tag 'multiqc'
  container 'quay.io/biocontainers/multiqc:1.31--pyhdfd78af_0'
  publishDir { "${params.outdir}/multiqc" }, mode: 'copy', overwrite: true
  stageInMode 'symlink'

  input:
  // Stage all QC files here as a barrier, then scan '.'
  path qc_files

  output:
  path 'multiqc_report.html', emit: report
  path 'multiqc_data',        emit: data
  path 'multiqc.log',         emit: log

  script:
  """
  set -euo pipefail
  {
    multiqc --version
    multiqc ${params.multiqc_extra_args ?: ''} -o . .
  } 2>&1 | tee multiqc.log

  test -s multiqc_report.html
  test -s multiqc.log
  """
}

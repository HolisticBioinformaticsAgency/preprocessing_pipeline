process HSMETRICS {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'
  publishDir { "${params.outdir}/${subject}/hsmetrics" }, mode: 'copy'

  input:
  tuple val(subject), val(sample_id), path(bam), path(bai), path(bed), path(ref_fa), path(ref_fai)

  output:
  tuple val(subject), val(sample_id), path("${sample_id}_hs_metrics.txt"), emit: hs

  script:
  def dict = "${ref_fa.baseName}.dict"
  """
  set -euo pipefail
  if [ ! -s "${dict}" ]; then
    picard CreateSequenceDictionary -R ${ref_fa} -O ${dict}
  fi
  picard BedToIntervalList -I ${bed} -O targets.interval_list -SD ${dict}
  picard CollectHsMetrics \
    -I ${bam} -O ${sample_id}_hs_metrics.txt -R ${ref_fa} \
    -BAIT_INTERVALS targets.interval_list \
    -TARGET_INTERVALS targets.interval_list \
    -VALIDATION_STRINGENCY SILENT
  test -s ${sample_id}_hs_metrics.txt
  """
}

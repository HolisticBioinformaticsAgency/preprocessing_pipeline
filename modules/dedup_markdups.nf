process DEDUP_MARKDUPS {
  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/picard:3.4.0--hdfd78af_0'

  // publish per subject
  publishDir { "${params.outdir}/${subject}/bams" }, mode: 'copy'

  input:
  // match upstream: (subject, sample_id, bam, bai)
  tuple val(subject), val(sample_id), path(bam), path(bai)

  output:
  // keep tuple shape; filenames reflect "marked duplicates"
  tuple val(subject),
        val(sample_id),
        path("${sample_id}_sorted_mkdup.bam"),
        path("${sample_id}_sorted_mkdup.bam.bai"),
        emit: bam

  // metrics for MultiQC etc.
  path "${sample_id}_dup_metrics.txt", emit: metrics

  script:
  """
  set -euo pipefail

  picard MarkDuplicates \
    I=${bam} \
    O=${sample_id}_sorted_mkdup.bam \
    M=${sample_id}_dup_metrics.txt \
    REMOVE_DUPLICATES=false \
    ASSUME_SORTED=true \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=SILENT

  # Picard may create ${sample_id}_sorted_mkdup.bai beside the BAM; normalize to .bam.bai if needed
  if [ -s "${sample_id}_sorted_mkdup.bai" ] && [ ! -e "${sample_id}_sorted_mkdup.bam.bai" ]; then
    ln -sf "${sample_id}_sorted_mkdup.bai" "${sample_id}_sorted_mkdup.bam.bai"
  fi

  # Sanity checks
  test -s "${sample_id}_sorted_mkdup.bam"
  test -s "${sample_id}_sorted_mkdup.bam.bai"
  test -s "${sample_id}_dup_metrics.txt"
  """
}
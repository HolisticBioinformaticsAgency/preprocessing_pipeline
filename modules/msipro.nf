/*
 * MSIPRO_PAIRED
 * Paired tumour/normal using `msisensor-pro msi`
 *
 * Input:
 *  (pub_base, subject, case_id,
 *   t_id, t_bam, t_bai,
 *   n_id, n_bam, n_bai,
 *   sites_tsv, fasta)
 *
 * Output:
 *  (pub_base, subject, case_id, <case_id>.msi)
 */
process MSIPRO_PAIRED {
  tag "${subject} (${case_id})"
  container 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'
  publishDir { "${pub_base}/msisensor" }, mode: 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id),
        val(t_id), path(tumor_bam),  path(tumor_bai),
        val(n_id), path(normal_bam), path(normal_bai),
        path(sites_tsv),
        path(fasta)

  output:
  tuple val(pub_base), val(subject), val(case_id), path("${case_id}.msi"), emit: msi

  script:
  """
  set -euo pipefail

  # Ensure BAM indexes exist alongside BAMs
  [ -e "${tumor_bam}.bai"  ] || ln -sf "${tumor_bai}"  "${tumor_bam}.bai"  || true
  [ -e "${normal_bam}.bai" ] || ln -sf "${normal_bai}" "${normal_bam}.bai" || true

  msisensor-pro msi \\
    -d ${sites_tsv} \\
    -n ${normal_bam} \\
    -t ${tumor_bam} \\
    -g ${fasta} \\
    -o ${case_id} \\
    -b ${task.cpus} \\
    > ${case_id}.log 2>&1

  # Normalize summary filename to <case_id>.msi
  if   [ -s "${case_id}.msi" ]; then
    :
  elif [ -s "${case_id}.msi.txt" ]; then
    mv "${case_id}.msi.txt" "${case_id}.msi"
  elif [ -s "${case_id}" ] && grep -q '^Total_Number_of_Sites' "${case_id}"; then
    mv "${case_id}" "${case_id}.msi"
  elif [ -s "${case_id}.txt" ] && grep -q '^Total_Number_of_Sites' "${case_id}.txt"; then
    mv "${case_id}.txt" "${case_id}.msi"
  else
    echo "[MSI] WARN: summary not found; emitting minimal placeholder" >&2
    echo -e "Total_Number_of_Sites\\t0\\nMSI_score\\tNA" > ${case_id}.msi
  fi

  test -s "${case_id}.msi"
  """
}

/*
 * MSIPRO_SINGLE
 * Single tumour (or single normal) using `msisensor-pro pro`
 *
 * Input:
 *  (pub_base, subject, case_id,
 *   t_id, t_bam, t_bai,
 *   sites_tsv, fasta)
 *
 * Output:
 *  (pub_base, subject, case_id, <case_id>.msi)
 */
process MSIPRO_SINGLE {
  tag "${subject} (${case_id})"
  container 'quay.io/biocontainers/msisensor-pro:1.3.0--hd979922_1'
  publishDir { "${pub_base}/msisensor" }, mode: 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id),
        val(t_id), path(tumor_bam), path(tumor_bai),
        path(sites_tsv),
        path(fasta)

  output:
  tuple val(pub_base), val(subject), val(case_id), path("${case_id}.msi"), emit: msi

  script:
  """
  set -euo pipefail

  # Ensure tumour BAM index exists alongside BAM
  [ -e "${tumor_bam}.bai" ] || ln -sf "${tumor_bai}" "${tumor_bam}.bai" || true

  msisensor-pro pro \\
    -d ${sites_tsv} \\
    -t ${tumor_bam} \\
    -g ${fasta} \\
    -o ${case_id} \\
    -b ${task.cpus} \\
    > ${case_id}.log 2>&1

  # Normalize summary filename to <case_id>.msi
  if   [ -s "${case_id}.msi" ]; then
    :
  elif [ -s "${case_id}.msi.txt" ]; then
    mv "${case_id}.msi.txt" "${case_id}.msi"
  elif [ -s "${case_id}" ] && grep -q '^Total_Number_of_Sites' "${case_id}"; then
    mv "${case_id}" "${case_id}.msi"
  elif [ -s "${case_id}.txt" ] && grep -q '^Total_Number_of_Sites' "${case_id}.txt"; then
    mv "${case_id}.txt" "${case_id}.msi"
  else
    echo "[MSI] WARN: summary not found; emitting minimal placeholder" >&2
    echo -e "Total_Number_of_Sites\\t0\\nMSI_score\\tNA" > ${case_id}.msi
  fi

  test -s "${case_id}.msi"
  """
}

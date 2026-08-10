/*
 * POLYSOLVER
 * Input:  (pub_base, subject, case_id, sample_id, bam, bai)
 * Output: (pub_base, subject, case_id, sample_id, <sample_id>_hla.txt)
 */
process POLYSOLVER {
  tag "${subject}_${sample_id}"
  container 'shub://IARCbioinfo/polysolver-singularity:v4'

  publishDir { "${pub_base}/polysolver" }, mode: 'copy'
  stageInMode 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id), val(sample_id), path(bam), path(bai)

  output:
  tuple val(pub_base), val(subject), val(case_id), val(sample_id), path("${sample_id}_hla.txt"), emit: hla

  script:
  def RACE       = params.polysolver_race       ?: 'Unknown'
  def BUILD      = params.polysolver_build      ?: 'hg38'
  def EMITVCF    = params.polysolver_emit_vcf   ?: 0
  def FQTYPE     = params.polysolver_fastqtype  ?: 'STDFQ'
  def INSERTCALC = params.polysolver_insertcalc ?: 0

  """
  set -euo pipefail

  WRAPPER="/home/polysolver/scripts/shell_call_hla_type"
  [ -x "\$WRAPPER" ] || { echo "ERROR: \$WRAPPER not found in container"; exit 127; }

  [ -e "${bam}.bai" ] || ln -sf "${bai}" "${bam}.bai" || true

  "\$WRAPPER" \\
    "${bam}" \\
    "${RACE}" \\
    "${EMITVCF}" \\
    "${BUILD}" \\
    "${FQTYPE}" \\
    "${INSERTCALC}" \\
    "\$(pwd)" > "${sample_id}_hla.txt"
  """
}

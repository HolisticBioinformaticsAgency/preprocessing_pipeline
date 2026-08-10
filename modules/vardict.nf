/*
 * VARDICT_SINGLE
 * Input:  (pub_base, subject, sample_id, bam, bai, ref_fa, ref_fai, bed)
 * Output:
 *   - (pub_base, subject, sample_id, <sample_id>_vardict.vcf)             [emit: vcf]
 */
process VARDICT_SINGLE {

  tag "${subject}_${sample_id}"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'

  input:
  tuple val(pub_base), val(subject), val(sample_id), path(bam), path(bai), path(ref_fa), path(ref_fai), path(bed)

  output:
  tuple val(pub_base), val(subject), val(sample_id), path("${sample_id}_vardict.vcf"), emit: vcf

  script:
  """
  set -euo pipefail
  export JAVA_TOOL_OPTIONS='-Xms2g -Xmx20g'

  # Determine VarDict -g column from BED (4 if present, else 0)
  GCOL=\$(awk 'BEGIN{c=0} !/^#/ && NF{print (NF>=4?4:0); exit}' ${bed})
  echo "[VarDict SINGLE] Using -g \${GCOL}" >&2

  vardict-java \\
    -G ${ref_fa} \\
    -N ${sample_id} \\
    -b ${bam} \\
    -th ${task.cpus} \\
    --nosv \\
    -f ${params.min_af} -z -c 1 -S 2 -E 3 -g \${GCOL} ${bed} \\
  | teststrandbias.R \\
  | var2vcf_valid.pl -A -N ${sample_id} -E -f ${params.min_af} \\
    > ${sample_id}_vardict.vcf
  """
}

/*
 * VARDICT_PAIRED
 * Input:  (pub_base, subject, case_id, t_id, t_bam, t_bai, n_id, n_bam, n_bai, ref_fa, ref_fai, bed)
 * Output:
 *   - (pub_base, subject, <case_id>_vardict.vcf)                           [emit: vcf]
 */
process VARDICT_PAIRED {

  tag "${subject} (T:${t_id}|N:${n_id})"
  container 'quay.io/biocontainers/vardict-java:1.8.3--hdfd78af_0'

  input:
  tuple val(pub_base), val(subject), val(case_id),
        val(t_id), path(t_bam), path(t_bai),
        val(n_id), path(n_bam), path(n_bai),
        path(ref_fa), path(ref_fai),
        path(bed)

  output:
  tuple val(pub_base), val(subject), val(case_id), path("${case_id}_vardict.vcf"), emit: vcf

  script:
  """
  set -euo pipefail
  export JAVA_TOOL_OPTIONS='-Xms2g -Xmx20g'

  # Determine VarDict -g column from BED (4 if present, else 0)
  GCOL=\$(awk 'BEGIN{c=0} !/^#/ && NF{print (NF>=4?4:0); exit}' ${bed})
  echo "[VarDict PAIRED] Using -g \${GCOL}" >&2

  TN="${t_id}|${n_id}"

  vardict-java \\
    -G ${ref_fa} \\
    -N "\${TN}" \\
    -b "${t_bam}|${n_bam}" \\
    -th ${task.cpus} \\
    --nosv \\
    -f ${params.min_af} -z -c 1 -S 2 -E 3 -g \${GCOL} ${bed} \\
  | testsomatic.R \\
  | var2vcf_paired.pl -A -N "\${TN}" -f ${params.min_af} \\
    > ${case_id}_vardict.vcf
  """
}

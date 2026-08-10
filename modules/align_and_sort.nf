process ALIGN_AND_SORT {
  tag "${subject}_${sample_id}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
      'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

  input:
  // (subject, sample_id, [R1,R2], staged FASTA, ABSOLUTE original FASTA path string)
  tuple val(subject), val(sample_id), path(reads), path(ref_fa), val(ref_src_abs)

  output:
  tuple val(subject), val(sample_id),
        path("${sample_id}_sorted.bam"),
        path("${sample_id}_sorted.bam.bai"),
        emit: bam

  script:
  """
  set -euo pipefail

  # Ensure BWA can find its sidecar index files next to the staged FASTA
  for ext in amb ann bwt pac sa; do
    if [ -s "${ref_src_abs}.\$ext" ]; then
      ln -sf "${ref_src_abs}.\$ext" .
    else
      echo "Missing BWA index: ${ref_src_abs}.\$ext" >&2
      exit 2
    fi
  done

  # Your RG and pipeline (no SAM written to disk)
  RG="@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPU:lib1\\tPL:Illumina"

  bwa mem -t ${task.cpus} -R "\$RG" ${ref_fa} ${reads.join(' ')} \\
    | samtools view -u -h - \\
    | samtools sort -@ ${task.cpus} -o "${sample_id}_sorted.bam"

  samtools index -@ ${task.cpus} "${sample_id}_sorted.bam" "${sample_id}_sorted.bam.bai"

  # Sanity checks
  test -s "${sample_id}_sorted.bam"
  test -s "${sample_id}_sorted.bam.bai"
  """
}
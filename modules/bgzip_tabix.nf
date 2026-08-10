process BGZIP_TABIX {
  tag "${subject}:${id}"
  container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
      'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/d7/d7e24dc1e4d93ca4d3a76a78d4c834a7be3985b0e1e56fddd61662e047863a8a/data' :
      'community.wave.seqera.io/library/bwa_htslib_samtools:83b50ff84ead50d0' }"

  publishDir { "${pub_base}/vcfs/${bucket}" }, mode: 'copy'

  input:
  tuple val(pub_base), val(bucket), val(subject), val(id), val(mode), path(vcf)

  output:
  tuple val(pub_base), val(bucket), val(subject), val(id), val(mode), path("${vcf.name}.gz"), emit: vcfgz
  path("${vcf.name}.gz.tbi"), emit: vcfgz_tbi

  script:
  """
  set -euo pipefail

  OUT_VCF="${vcf.name}.gz"
  bgzip -f -c "${vcf}" > "\${OUT_VCF}"
  tabix -f -p vcf "\${OUT_VCF}"
  """
}

process BCFTOOLS_NORM {
  tag "${subject}:${id}"
  publishDir { "${pub_base}/vcfs/norm" }, mode: 'copy'
  container 'quay.io/biocontainers/bcftools:1.23.1--hb2cee57_0'

  input:
  tuple val(pub_base), val(subject), val(id), val(mode), path(vcf), path(ref_fa), path(ref_fai)

  output:
  tuple val(pub_base), val(subject), val(id), val(mode),
        path("${vcf.simpleName}_norm.vcf.gz"),
        path("${vcf.simpleName}_norm.vcf.gz.tbi"), emit: vcfgz

  script:
  def out = "${vcf.simpleName}_norm.vcf.gz"
  """
  set -euo pipefail
  bcftools reheader -f ${ref_fai} ${vcf} \
  | bcftools norm \
      -f ${ref_fa} \
      -c wx \
      -m -any \
      -O v \
  | bgzip > ${out}
  tabix -p vcf ${out}
  """
}

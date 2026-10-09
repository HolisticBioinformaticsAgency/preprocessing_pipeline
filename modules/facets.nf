process FACETS {
  tag "${subject}_${case_id}"
  container 'quay.io/biocontainers/cnv_facets:0.16.1--py312r43h8537716_1'

  publishDir { "${pub_base}/facets_purity" }, mode: 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id),
        val(t_id), path(t_bam, stageAs: 'tumour.bam'), path(t_bai, stageAs: 'tumour.bam.bai'),
        val(n_id), path(n_bam, stageAs: 'normal.bam'), path(n_bai, stageAs: 'normal.bam.bai')

  output:
  tuple val(case_id), path("${case_id}.facets_purity.txt"),  emit: purity
  tuple val(case_id), path("${case_id}.facets_ploidy.txt"),  emit: ploidy
  path "${case_id}.facets.tsv",                              emit: summary
  path "${case_id}.snp_pileup.csv.gz",                       emit: pileup
  path "${case_id}*.cnv.pdf",    optional: true
  path "${case_id}*.spider.pdf", optional: true
  path "${case_id}*.fit.rds",    optional: true

  script:
  """
  set -euo pipefail

  SNPVCF='${params.common_snps}'
  test -s "\$SNPVCF"

  # 1. snp-pileup — normal BAM first, then tumour (Tu parameters)
  snp-pileup -g \\
    --min-map-quality  ${params.min_map_quality} \\
    --min-base-quality ${params.min_base_quality} \\
    --pseudo-snps      ${params.pseudo_snps} \\
    --max-depth        ${params.max_depth} \\
    "\$SNPVCF" ${case_id}.snp_pileup.csv.gz normal.bam tumour.bam

  # 2. auto FACETS fit  (facets_fit.R is staged from bin/)
  facets_fit.R ${case_id}.snp_pileup.csv.gz ${case_id} . hg38 auto

  PLOIDY=\$(awk -F'\\t' 'NR==2{print \$5}' ${case_id}.facets.tsv)
  FLAG=\$(awk   -F'\\t' 'NR==2{print \$8}' ${case_id}.facets.tsv)
  echo "[FACETS] ${case_id} auto fit: ploidy=\${PLOIDY} flag=\${FLAG}"

  # 3. always use the auto fit — no dipLogR=0 refit
  CHOSEN=${case_id}.facets.tsv

  # 4. chosen purity (column 4) -> single-value file for CNVKIT (--purity)
  awk -F'\\t' 'NR==2 {print (\$4 == "" ? "NA" : \$4)}' "\${CHOSEN}" > ${case_id}.facets_purity.txt
  echo "[FACETS] ${case_id} chosen purity = \$(cat ${case_id}.facets_purity.txt)  (from \${CHOSEN})"

  # 5. chosen ploidy (column 5) -> single-value file for CNVKIT (--ploidy)
  awk -F'\\t' 'NR==2 {print (\$5 == "" ? "NA" : \$5)}' "\${CHOSEN}" > ${case_id}.facets_ploidy.txt
  echo "[FACETS] ${case_id} chosen ploidy = \$(cat ${case_id}.facets_ploidy.txt)  (from \${CHOSEN})"
  """
}

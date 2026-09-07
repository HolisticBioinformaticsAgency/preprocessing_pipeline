process CNVKIT {
  tag "${subject}_${case_id}"
  container 'quay.io/biocontainers/cnvkit:0.9.12--pyhdfd78af_1'

  publishDir { "${pub_base}/cnvkit" }, mode: 'copy'

  input:
  tuple val(pub_base), val(subject), val(case_id), val(mode),
        val(main_id),
        path(main_bam, stageAs: 'case.bam'),
        path(main_bai, stageAs: 'case.bam.bai'),
        val(ctrl_id),
        path(ctrl_bam, stageAs: 'control.bam'),
        path(ctrl_bai, stageAs: 'control.bam.bai'),
        path(norm_vcf, stageAs: 'norm.vcf.gz'),
        path(norm_tbi, stageAs: 'norm.vcf.gz.tbi')

  output:
  tuple val(pub_base), val(subject), val(case_id),
        path("${case_id}.cnr"),
        path("${case_id}.cbs.dlc.cns"),
        emit: cnv
  path "${case_id}.cbs.dlc.call.cns",    emit: call_cns
  path "${case_id}.genemetrics_cnr.tsv", emit: genemetrics_cnr
  path "${case_id}.genemetrics_cns.tsv", emit: genemetrics_cns
  path "${case_id}-scatter.pdf", optional: true, emit: scatter
  path "${case_id}-diagram.pdf", optional: true, emit: diagram

  script:
  def call_cmd = (mode == 'paired')
    ? """
  cnvkit.py call ${case_id}.cbs.dlc.cns \\
    --vcf ${norm_vcf} \\
    --sample-id ${main_id} \\
    --normal-id ${ctrl_id} \\
    --zygosity-freq 0.35 \\
    --min-variant-depth 30 \\
    --output ${case_id}.cbs.dlc.call.cns
"""
    : """
  cnvkit.py call ${case_id}.cbs.dlc.cns \\
    --output ${case_id}.cbs.dlc.call.cns
"""
  def scatter_vcf = (mode == 'paired')
    ? "--vcf ${norm_vcf} --sample-id ${main_id} --normal-id ${ctrl_id} --zygosity-freq 0.35 --min-variant-depth 30"
    : ""
  """
  set -euo pipefail
  export MPLCONFIGDIR='/tmp'
  export MPLBACKEND='Agg'
  export XDG_CACHE_HOME='/tmp'
  export HOME='/tmp'

  # Build .cnr from BAM against the panel of normals
  cnvkit.py batch case.bam \\
    --method hybrid \\
    --segment-method none \\
    --reference ${params.cnvkit_pon} \\
    -p ${task.cpus} \\
    -d .

  mv case.cnr ${case_id}.cnr

  # CBS segmentation
  cnvkit.py segment ${case_id}.cnr \\
    --method cbs \\
    --drop-low-coverage \\
    -o ${case_id}.cbs.dlc.cns

  ${call_cmd}

  # Gene-level metrics: on the .cnr alone, and on the .cnr + called .cns
  cnvkit.py genemetrics ${case_id}.cnr \\
    --drop-low-coverage \\
    -t 0 \\
    -o ${case_id}.genemetrics_cnr.tsv

  cnvkit.py genemetrics ${case_id}.cnr \\
    -s ${case_id}.cbs.dlc.call.cns \\
    -t 0 \\
    -o ${case_id}.genemetrics_cns.tsv

  # Whole-genome scatter + diagram (use the called .cns so segments are
  # coloured by CNA status / LOH; scatter also gets the BAF panel when paired)
  cnvkit.py scatter ${case_id}.cnr -s ${case_id}.cbs.dlc.call.cns ${scatter_vcf} -o ${case_id}-scatter.pdf || true
  cnvkit.py diagram ${case_id}.cnr -s ${case_id}.cbs.dlc.call.cns -o ${case_id}-diagram.pdf || true
  """
}

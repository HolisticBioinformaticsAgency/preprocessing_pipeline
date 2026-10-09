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
        path(norm_tbi, stageAs: 'norm.vcf.gz.tbi'),
        val(facets_purity),
        val(facets_ploidy),
        val(sample_sex)

  output:
  tuple val(pub_base), val(subject), val(case_id),
        path("${case_id}.cnr"),
        path("${case_id}.cbs.dlc.cns"),
        emit: cnv
  path "${case_id}.cbs.dlc.call.cns",              emit: call_cns
  path "${case_id}.cbs.dlc.purity.ploidy.call.cns", optional: true, emit: purity_call_cns
  path "${case_id}.genemetrics_cnr.tsv",           emit: genemetrics_cnr
  path "${case_id}.genemetrics_cns.tsv",           emit: genemetrics_cns
  path "${case_id}.genemetrics_call_cns.tsv",      emit: genemetrics_call_cns
  path "${case_id}-scatter.pdf", optional: true,   emit: scatter
  path "${case_id}-diagram.pdf", optional: true,   emit: diagram

  script:
  // BAF/allele-specific args shared by the paired call, purity call and scatter
  def paired_vcf_args = "--vcf ${norm_vcf} --sample-id ${main_id} --normal-id ${ctrl_id} --zygosity-freq 0.35 --min-variant-depth 30"

  // Chromosomal sex from the samplesheet, translated to CNVkit's own -x/--sample-sex
  // vocabulary (lowercase 'm'/'f' — 'M'/'F'/'NA'/blank do not match its argparse
  // choices). Left blank ('') when unknown, so CNVkit falls back to its own
  // coverage-based auto-detection, unchanged from today's behaviour.
  def sex_norm = (sample_sex ?: '').toString().trim().toUpperCase()
  def sex_arg  = (sex_norm == 'M') ? '-x male' : (sex_norm == 'F') ? '-x female' : ''

  // Plain call (no purity) — always produced
  def call_cmd = (mode == 'paired')
    ? """
  cnvkit.py call ${case_id}.cbs.dlc.cns \\
    ${paired_vcf_args} \\
    ${sex_arg} \\
    --output ${case_id}.cbs.dlc.call.cns
"""
    : """
  cnvkit.py call ${case_id}.cbs.dlc.cns \\
    ${sex_arg} \\
    --output ${case_id}.cbs.dlc.call.cns
"""

  // Purity call — only when FACETS returned a numeric purity (paired cases)
  def has_purity = (mode == 'paired' && facets_purity && facets_purity != 'NA' && (facets_purity ==~ /^[0-9]*\.?[0-9]+$/))

  // Ploidy from FACETS, rounded to the nearest integer (CNVkit's --ploidy is int-typed).
  // Only applied alongside a valid purity, on paired cases — ploidy correction without
  // purity correction isn't a combination CNVkit's calling model is used for here.
  def has_ploidy  = (has_purity && facets_ploidy && facets_ploidy != 'NA' && (facets_ploidy ==~ /^[0-9]*\.?[0-9]+$/))
  def ploidy_arg  = has_ploidy ? "--ploidy ${Math.round(facets_ploidy.toDouble())}" : ''

  def purity_call_cmd = has_purity
    ? """
  cnvkit.py call ${case_id}.cbs.dlc.cns \\
    ${paired_vcf_args} \\
    --purity ${facets_purity} \\
    ${ploidy_arg} \\
    ${sex_arg} \\
    --output ${case_id}.cbs.dlc.purity.ploidy.call.cns
"""
    : """
  echo "[CNVKIT] ${case_id}: no FACETS purity (${facets_purity}) - skipping --purity call"
"""

  // Downstream (genemetrics -s, scatter, diagram) use the purity call when it exists
  def eff_cns = has_purity ? "${case_id}.cbs.dlc.purity.ploidy.call.cns" : "${case_id}.cbs.dlc.call.cns"

  def scatter_vcf = (mode == 'paired') ? paired_vcf_args : ""
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
  ${purity_call_cmd}

  # Gene-level metrics: on the .cnr alone, and on the .cnr + the effective called .cns
  cnvkit.py genemetrics ${case_id}.cnr \\
    --drop-low-coverage \\
    -t 0 \\
    -o ${case_id}.genemetrics_cnr.tsv

  cnvkit.py genemetrics ${case_id}.cnr \\
    -s ${eff_cns} \\
    -t 0 \\
    -o ${case_id}.genemetrics_cns.tsv

  # The same on the plain (no purity) call: CGenGx reads the plain call, so its gene table has to come from it too (CGen, 2026-10-09)
  cnvkit.py genemetrics ${case_id}.cnr \\
    -s ${case_id}.cbs.dlc.call.cns \\
    -t 0 \\
    -o ${case_id}.genemetrics_call_cns.tsv

  # Whole-genome scatter + diagram (use the effective called .cns so segments are
  # coloured by CNA status / LOH; scatter also gets the BAF panel when paired)
  cnvkit.py scatter ${case_id}.cnr -s ${eff_cns} ${scatter_vcf} -o ${case_id}-scatter.pdf || true
  cnvkit.py diagram ${case_id}.cnr -s ${eff_cns} -o ${case_id}-diagram.pdf || true
  """
}

workflow INIT_PARAMS {
  if( !params.containsKey('samplesheet')   ) params.samplesheet    = "./refs/samplesheet.csv"
  if( !params.containsKey('reference')     ) params.reference      = null
  if( !params.containsKey('reference_fai') ) params.reference_fai  = null
  if( !params.containsKey('bed')           ) params.bed            = null
  if( !params.containsKey('outdir')        ) params.outdir         = "results"
  if( !params.containsKey('min_af')        ) params.min_af         = 0.05
  if( !params.containsKey('multiqc_extra_args') ) params.multiqc_extra_args = ''

  // MSI tumor-only baseline (optional)
  if( !params.containsKey('msipro_baseline') ) params.msipro_baseline = null

  // Polysolver defaults
  if( !params.containsKey('polysolver_race')       ) params.polysolver_race       = 'Caucasian'
  if( !params.containsKey('polysolver_build')      ) params.polysolver_build      = 'hg38'
  if( !params.containsKey('polysolver_emit_vcf')   ) params.polysolver_emit_vcf   = 0
  if( !params.containsKey('polysolver_fastqtype')  ) params.polysolver_fastqtype  = 'STDFQ'
  if( !params.containsKey('polysolver_insertcalc') ) params.polysolver_insertcalc = 0
}

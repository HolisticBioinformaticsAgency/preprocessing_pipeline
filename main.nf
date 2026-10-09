// -------------------- Includes --------------------
include { INIT_PARAMS }                    from './modules/init_params.nf'
include { FASTQC }                         from './modules/fastqc.nf'
include { ALIGN_AND_SORT }                 from './modules/align_and_sort.nf'
include { DEDUP_MARKDUPS }                 from './modules/dedup_markdups.nf'
include { HSMETRICS }                      from './modules/hsmetrics.nf'
include { POLYSOLVER }                     from './modules/polysolver.nf'

include { VARDICT_SINGLE }                 from './modules/vardict.nf'
include { VARDICT_PAIRED }                 from './modules/vardict.nf'

include { FACETS }                         from './modules/facets.nf'
include { CNVKIT }                         from './modules/cnvkit.nf'

include { MSIPRO_SCAN }                    from './modules/msipro_scan.nf'
include { MSIPRO_PAIRED ; MSIPRO_SINGLE }  from './modules/msipro.nf'

include { MULTIQC }                        from './modules/multiqc.nf'
include { BCFTOOLS_NORM }                  from './modules/bcftools_norm.nf'

// ---- same module, multiple aliases ----
include { BGZIP_TABIX as BGZIP_VARDICT_SINGLE } from './modules/bgzip_tabix.nf'
include { BGZIP_TABIX as BGZIP_VARDICT_PAIRED } from './modules/bgzip_tabix.nf'

// -------------------- Workflow --------------------
workflow {

  // -------------------- Samplesheet → Channels --------------------
  def REQUIRED = ['sample','subject','status','fastq_1','fastq_2']

  Channel
    .fromPath(params.samplesheet)
    .splitCsv(header:true)
    .map { row ->
      def missing = REQUIRED.findAll { k ->
        !row.containsKey(k) || row[k] == null || row[k].toString().trim() == ''
      }
      if( missing ) {
        throw new IllegalArgumentException("Samplesheet row has missing/empty field(s): ${missing} → ${row}")
      }
      tuple(
        row.sample as String,
        row.subject as String,
        (row.status as String).toLowerCase(),
        file(row.fastq_1),
        file(row.fastq_2),
        ((row.sex ?: 'NA') as String)
      )
    }
    .set { ch_sheet }

  // Reads for QC/BWA: (subject, sample, [R1,R2])
  ch_reads = ch_sheet.map { sample, subject, status, r1, r2, sex ->
    tuple(subject, sample, [ r1, r2 ])
  }

  // Sample metadata for joins: (sample, subject, status, sex)
  ch_meta  = ch_sheet.map { sample, subject, status, r1, r2, sex ->
    tuple(sample, subject, status, sex)
  }

  // Reference & BED singletons
  Channel.of( file(params.reference)     ).set { ch_ref_fa }
  Channel.of( file(params.reference_fai) ).set { ch_ref_fai }
  Channel.of( file(params.bed)           ).set { ch_bed }

  // Also pass original absolute FASTA path (for BWA sidecars)
  Channel
    .of( params.reference )
    .map { new File(it as String).getAbsolutePath() }
    .set { ch_ref_src_abs }

  // ---------- QC ----------
  INIT_PARAMS()
  FASTQC( ch_reads )

  // ---------- Align + Sort ----------
  ch_align_in   = ch_reads.combine(ch_ref_fa).combine(ch_ref_src_abs)
  ch_bam_sorted = ALIGN_AND_SORT( ch_align_in )

  // ---------- Dedup ----------
  (ch_bam, ch_dedup_metrics) = DEDUP_MARKDUPS( ch_bam_sorted )

  // ---------- Attach metadata ----------
  ch_bam_meta = ch_bam
    .map  { sub, sample, bam, bai -> tuple(sample, tuple(sub, bam, bai)) }
    .join ( ch_meta )
    .map  { sample, left, subject, status, sex ->
      tuple(subject, sample, status, sex, left[1], left[2])
    }

  // ---------- HsMetrics ----------
  def ch_hs_in = ch_bam.combine(ch_bed).combine(ch_ref_fa).combine(ch_ref_fai)
    .map { sub, sample, bam, bai, bed, ref_fa, ref_fai -> tuple(sub, sample, bam, bai, bed, ref_fa, ref_fai) }
  HSMETRICS( ch_hs_in )

  // ---------- Group by subject ----------
  def by_subject = ch_bam_meta
    .map { sub, sample, status, sex, bam, bai -> tuple(sub, tuple(sample, status, bam, bai, sex)) }
    .groupTuple()
    .map { sub, recs -> tuple(sub, recs.sort { a, b -> a[0] <=> b[0] }) }

  // ---------- Compute cases and publishing base ----------
  // Each sample record carries its samplesheet 'sex' as a trailing element (index 4),
  // appended after the existing (sample, bam, bai, role) shape so nothing downstream
  // that only reads indices 0-3 needs to change.
  def ch_cases_pub = by_subject.flatMap { sub, recs ->
    def tumors  = recs.findAll { it[1] == 'tumor'  }.collect { [ it[0], it[2], it[3], 'tumor',  it[4] ] }
    def normals = recs.findAll { it[1] == 'normal' }.collect { [ it[0], it[2], it[3], 'normal', it[4] ] }

    def ncases = (tumors && normals) ? (tumors.size() * normals.size()) : recs.size()

    def out = []
    if( tumors && normals ) {
      tumors.each { t ->
        normals.each { n ->
          def case_id = "${t[0]}_${n[0]}"
          def base    = (ncases > 1) ? "${params.outdir}/${sub}/${case_id}" : "${params.outdir}/${sub}"
          out << tuple(sub, case_id, 'paired', [t, n], base)
        }
      }
    } else {
      recs.each { r ->
        def case_id = r[0]
        def base    = (ncases > 1) ? "${params.outdir}/${sub}/${case_id}" : "${params.outdir}/${sub}"
        out << tuple(sub, case_id, 'single', [ [ r[0], r[2], r[3], r[1], r[4] ] ], base)
      }
    }
    out
  }

  // ---------- POLYSOLVER ----------
  POLYSOLVER(
    ch_cases_pub.flatMap { sub, case_id, mode, samples, pub_base ->
      samples.collect { s ->
        tuple(pub_base, sub, case_id, s[0], s[1], s[2])
      }
    }
  )

  // ---------- VarDict in ----------
  def ch_vardict_single_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'single' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .map { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def s = samples[0]
      tuple(pub_base, sub, s[0], s[1], s[2], ref_fa, ref_fai, bed)
    }

  def ch_vardict_paired_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine(ch_ref_fa).combine(ch_ref_fai).combine(ch_bed)
    .flatMap { sub, case_id, mode, samples, pub_base, ref_fa, ref_fai, bed ->
      def tumors  = samples.findAll { it[3] == 'tumor'  }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], ref_fa, ref_fai, bed)
        }
      }
    }

  // ---------- VarDict run ----------
  VARDICT_SINGLE( ch_vardict_single_in )
  def ch_vardict_single_vcf_raw = VARDICT_SINGLE.out.vcf
    .map { pub, sub, id, vcf -> tuple(pub, 'vardict', sub, id, 'single', vcf) }

  VARDICT_PAIRED( ch_vardict_paired_in )
  def ch_vardict_paired_vcf_raw = VARDICT_PAIRED.out.vcf
    .map { pub, sub, case_id, vcf -> tuple(pub, 'vardict', sub, case_id, 'paired', vcf) }

  // ---------- bgzip/tabix VarDict ----------
  BGZIP_VARDICT_SINGLE( ch_vardict_single_vcf_raw )
  BGZIP_VARDICT_PAIRED( ch_vardict_paired_vcf_raw )

  // ---------- bcftools norm ----------
  def ch_vardict_single_vcfgz = BGZIP_VARDICT_SINGLE.out.vcfgz
    .map { pub, bucket, sub, id, mode, vcf -> tuple(pub, sub, id, mode, vcf) }

  def ch_vardict_paired_vcfgz = BGZIP_VARDICT_PAIRED.out.vcfgz
    .map { pub, bucket, sub, id, mode, vcf -> tuple(pub, sub, id, mode, vcf) }

  def ch_norm_in = ch_vardict_single_vcfgz.mix( ch_vardict_paired_vcfgz )
    .combine( ch_ref_fa )
    .combine( ch_ref_fai )
  BCFTOOLS_NORM( ch_norm_in )

  // ---------- FACETS (paired cases only) ----------
  // Tumour purity/ploidy from the matched pair; feeds CNVkit `call --purity`.
  def ch_facets_in = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .flatMap { sub, case_id, mode, samples, pub_base ->
      def tumors  = samples.findAll { it[3] == 'tumor' }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2])
        }
      }
    }

  FACETS( ch_facets_in )

  // (case_id, purity-string)   purity-string is a number or "NA"
  def ch_facets_purity = FACETS.out.purity
    .map { case_id, pfile -> tuple(case_id, pfile.text.trim()) }

  // (case_id, ploidy-string)   ploidy-string is a number or "NA"
  def ch_facets_ploidy = FACETS.out.ploidy
    .map { case_id, pfile -> tuple(case_id, pfile.text.trim()) }

  // ---------- CNVkit per case ----------
  // Norm'd VarDict VCF for each case, keyed by case_id, used by CNVkit call for BAF/allele-specific CN
  def ch_norm_vcf = BCFTOOLS_NORM.out.vcfgz
    .map { pub, sub, id, mode, vcfgz, tbi -> tuple(id, vcfgz, tbi) }

  def ch_cnv_jobs = ch_cases_pub
    .flatMap { sub, case_id, mode, samples, pub_base ->
      if( mode == 'paired' ) {
        def tumors  = samples.findAll { it[3] == 'tumor' }
        def normals = samples.findAll { it[3] == 'normal' }
        tumors.collectMany { t ->
          normals.collect { n ->
            tuple(case_id, pub_base, sub, mode, t[0], t[1], t[2], n[0], n[1], n[2], t[4])
          }
        }
      } else {
        def s = samples[0]
        [ tuple(case_id, pub_base, sub, mode,
                s[0], s[1], s[2],
                'NA', file('/dev/null'), file('/dev/null'), s[4]) ]
      }
    }
    .join( ch_norm_vcf )
    // left-join FACETS purity/ploidy: single-mode / FACETS-less cases get null -> 'NA'
    .join( ch_facets_purity, remainder: true )
    .join( ch_facets_ploidy, remainder: true )
    .map { case_id, pub_base, sub, mode, main_id, main_bam, main_bai, ctrl_id, ctrl_bam, ctrl_bai, sex, vcfgz, tbi, purity, ploidy ->
      tuple(pub_base, sub, case_id, mode, main_id, main_bam, main_bai, ctrl_id, ctrl_bam, ctrl_bai, vcfgz, tbi, (purity ?: 'NA'), (ploidy ?: 'NA'), (sex ?: 'NA'))
    }

  CNVKIT( ch_cnv_jobs )
  def ch_cnv     = CNVKIT.out.cnv
  def ch_scatter = CNVKIT.out.scatter
  def ch_diagram = CNVKIT.out.diagram // do we need these here?

  // ---------- MSI ----------
  MSIPRO_SCAN( ch_ref_fa )
  def ch_msisites  = MSIPRO_SCAN.out.sites
  def ch_msiref_fa = MSIPRO_SCAN.out.fasta

  def ch_msi_jobs_paired = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'paired' }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .flatMap { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def tumors  = samples.findAll { it[3] == 'tumor' }
      def normals = samples.findAll { it[3] == 'normal' }
      tumors.collectMany { t ->
        normals.collect { n ->
          tuple(pub_base, sub, case_id, t[0], t[1], t[2], n[0], n[1], n[2], sites_tsv, ref_fa)
        }
      }
    }

  def ch_msi_jobs_single_tumor = ch_cases_pub
    .filter { sub, case_id, mode, samples, pub_base -> mode == 'single' && samples[0][3] == 'tumor' }
    .combine( ch_msisites )
    .combine( ch_msiref_fa )
    .map { sub, case_id, mode, samples, pub_base, sites_tsv, ref_fa ->
      def t = samples[0]
      tuple(pub_base, sub, case_id, t[0], t[1], t[2], sites_tsv, ref_fa)
    }

  MSIPRO_PAIRED ( ch_msi_jobs_paired )
  MSIPRO_SINGLE ( ch_msi_jobs_single_tumor )

  // ---------- MultiQC ----------
  def ch_hsmetrics_files = HSMETRICS.out.hs.map { sub, sid, hs -> hs }
  def ch_multiqc_inputs = ch_dedup_metrics.mix(ch_hsmetrics_files).mix(FASTQC.out.zip).collect()
  MULTIQC( ch_multiqc_inputs )
}


// Run summary for the copy-number record and any report built on these outputs (CGen, 2026-10-09):
// code version, Nextflow version, dates, outcome, command line and every parameter.
workflow.onComplete {
    def rev = workflow.commitId ?: "git -C ${projectDir} rev-parse HEAD".execute().text.trim()
    def dirty = workflow.commitId ? "" : "git -C ${projectDir} status --porcelain".execute().text.trim()
    def info = [revision: rev, uncommitted_changes: dirty ? true : false, nextflow: nextflow.version.toString(),
                started: workflow.start.toString(), completed: workflow.complete.toString(), success: workflow.success,
                command: workflow.commandLine, params: params]
    file("${params.outdir}/pipeline_info").mkdirs()
    file("${params.outdir}/pipeline_info/run_info_${workflow.start.format('yyyyMMdd_HHmmss')}.json").text =
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(info))
}

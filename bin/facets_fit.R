#!/usr/bin/env Rscript
suppressMessages(library(facets))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5L)
  stop("usage: facets_fit.R <pileup> <sample> <outdir> <gbuild> <auto|number>")

pileup <- args[1]
sample <- args[2]
outdir <- args[3]
gbuild <- args[4]
dparg  <- args[5]

## --- facets parameters --------------------------------------------------------
## facets defaults; only the snp-pileup parameters are changed for this re-run
## (those are set in the shell driver). Kept here as named constants so a future
## tweak is a one-line edit.
SEED      <- 1234L
NDEPTH    <- 35L     # preProcSample: min normal depth to keep a position
NDEPTHMAX <- 1000L   # preProcSample: max normal depth to keep a position
SNP_NBHD  <- 250L    # preProcSample: SNP-thinning neighbourhood
CVAL_PRE  <- 25L     # preProcSample: pre-segmentation critical value
CVAL_PROC <- 150L    # procSample: segmentation critical value (exome)
MIN_NHET  <- 15L     # procSample / emcncf: min het SNPs per segment

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(SEED)

if (identical(dparg, "auto")) {
  dipLogR <- NULL
  tag     <- ""
  mode    <- "auto"
} else {
  dipLogR <- suppressWarnings(as.numeric(dparg))
  if (is.na(dipLogR)) stop("<diplogr> must be 'auto' or a number")
  tag  <- paste0(".diplogr", dparg)
  mode <- paste0("fixed:", dparg)
}

out_tsv    <- file.path(outdir, paste0(sample, tag, ".facets.tsv"))
out_rds    <- file.path(outdir, paste0(sample, tag, ".fit.rds"))
out_pdf    <- file.path(outdir, paste0(sample, tag, ".cnv.pdf"))
out_spider <- file.path(outdir, paste0(sample, tag, ".spider.pdf"))

write_row <- function(purity, ploidy, dlr, nhet, nseg, flag) {
  num <- function(v, d = 4L) if (is.null(v) || length(v) == 0L || is.na(v)) NA else round(as.numeric(v), d)
  df <- data.frame(
    sample        = sample,
    diplogr_mode  = mode,
    diplogr_value = num(dlr),
    purity        = num(purity),
    ploidy        = num(ploidy),
    nhet          = if (is.null(nhet) || is.na(nhet)) NA else as.integer(nhet),
    nseg          = if (is.null(nseg) || is.na(nseg)) NA else as.integer(nseg),
    flag          = flag,
    stringsAsFactors = FALSE
  )
  write.table(df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
}

invisible(tryCatch({
  rcmat <- readSnpMatrix(pileup)
  xx <- preProcSample(rcmat, gbuild = gbuild,
                      ndepth = NDEPTH, ndepthmax = NDEPTHMAX,
                      snp.nbhd = SNP_NBHD, cval = CVAL_PRE)
  oo  <- procSample(xx, cval = CVAL_PROC, min.nhet = MIN_NHET, dipLogR = dipLogR)
  fit <- emcncf(oo, min.nhet = MIN_NHET)

  nseg <- nrow(oo$out)
  nhet <- sum(oo$out$nhet, na.rm = TRUE)

  try({ pdf(out_pdf, width = 9, height = 7); plotSample(x = oo, emfit = fit); dev.off() }, silent = TRUE)
  try({ pdf(out_spider, width = 6, height = 6); logRlogORspider(oo$out, oo$dipLogR); dev.off() }, silent = TRUE)
  saveRDS(list(oo = oo, fit = fit, params = list(
    ndepth = NDEPTH, ndepthmax = NDEPTHMAX, snp.nbhd = SNP_NBHD,
    cval.pre = CVAL_PRE, cval.proc = CVAL_PROC, min.nhet = MIN_NHET,
    gbuild = gbuild, diplogr_mode = mode)), out_rds)

  write_row(fit$purity, fit$ploidy, oo$dipLogR, nhet, nseg, "ok")
  cat(sprintf("[facets_fit] %s (%s): purity=%s ploidy=%s dipLogR=%s -> %s\n",
              sample, mode,
              ifelse(is.null(fit$purity) || is.na(fit$purity), "NA", round(fit$purity, 3)),
              ifelse(is.null(fit$ploidy) || is.na(fit$ploidy), "NA", round(fit$ploidy, 3)),
              round(oo$dipLogR, 3), out_tsv))
}, error = function(e) {
  write_row(NA, NA, dipLogR, NA, NA, paste0("ERROR: ", conditionMessage(e)))
  cat(sprintf("[facets_fit] %s (%s): ERROR: %s\n", sample, mode, conditionMessage(e)))
}))

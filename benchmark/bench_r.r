#!/usr/bin/env Rscript
# bench_r.R
# Benchmarks R GWAS visualization tools: CMplot, qqman.
#
# Usage:
#   Rscript bench_r.R --tool CMplot --input data/sumstats_1M.tsv \
#       --size 1M --replicates 5 --outdir results/ --figdir figures/
#
# Writes one CSV row per replicate to results/bench_r.csv.

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
})

# ---------------------------------------------------------------------------
# CLI arguments
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--tool",       type="character", help="Tool: CMplot or qqman"),
  make_option("--input",      type="character", help="Path to sumstats TSV"),
  make_option("--size",       type="character", help="Dataset size label (e.g. 1M)"),
  make_option("--plot-type",  type="character", default="manhattan",
              help="Plot type: manhattan, circular, qq [default: manhattan]"),
  make_option("--replicates", type="integer",   default=5L),
  make_option("--manifest",   type="character", default=NULL,
              help="Path to manifest file listing trait TSV paths, one per line.\n\t\t    Used for multitrack_manhattan and multitrack_circular plot types.\n\t\t    The column-wise merge required by CMplot is performed inside\n\t\t    run_CMplot_multitrack() as part of the timed section."),
  make_option("--replicate",  type="integer",   default=NULL,
              help="Run only this single replicate (1-based). Used by bench_r_timed.sh\n\t\t    to isolate each rep in its own process for /usr/bin/time -v memory capture.\n\t\t    When set, peak_mem_mb is written as EXTERNAL (filled in by the wrapper)."),
  make_option("--outdir",     type="character", default="results"),
  make_option("--figdir",     type="character", default="figures")
)

opt <- parse_args(OptionParser(option_list=option_list))
tool       <- opt$tool
input      <- opt$input
size       <- opt$size
plot_type  <- opt[["plot-type"]]
reps       <- opt$replicates
single_rep <- opt$replicate   # NULL unless --replicate N supplied
manifest   <- opt$manifest    # NULL unless --manifest path supplied
outdir     <- opt$outdir
figdir     <- opt$figdir

dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
dir.create(figdir, showWarnings=FALSE, recursive=TRUE)

# ---------------------------------------------------------------------------
# Tool wrappers
# ---------------------------------------------------------------------------

run_CMplot <- function(df, out_path, plot_type) {
  suppressPackageStartupMessages(library(CMplot))

  # CMplot expects: SNP, CHR, BP, then one or more p-value columns.
  # This specific column order is a known usability friction point.
  cm_df <- df[, c("SNP", "CHR", "BP", "P")]

  # --- CMplot output file naming fix ---------------------------------------
  # CMplot constructs its own filename by prepending a plot-type prefix
  # (e.g. "Rect_Manhtn.", "Circ_Manhtn.", "QQplot.") to whatever you pass
  # as file.name, and writes the result relative to the *working directory*.
  # Passing a full absolute path as file.name therefore produces an invalid
  # concatenation like "Rect_Manhtn./full/path/stem" which cannot be opened.
  #
  # Fix: (1) cd into the figures directory, (2) pass only the bare stem,
  # (3) use list.files() with a pattern to find whatever CMplot actually
  #     wrote (prefix varies across CMplot versions), then rename to out_path.
  # -------------------------------------------------------------------------
  out_dir  <- normalizePath(dirname(out_path), mustWork = TRUE)
  out_stem <- tools::file_path_sans_ext(basename(out_path))

  # Snapshot files present before calling CMplot so we can identify new ones
  before <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)

  # Temporarily change working directory; restore on function exit
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  if (plot_type == "manhattan") {
    CMplot(
      cm_df,
      plot.type   = "m",
      LOG10       = TRUE,
      threshold   = c(5e-8, 1e-5),
      file        = "png",
      file.name   = out_stem,       # bare stem only — no path, no extension
      file.output = TRUE,
      verbose     = FALSE,
      width       = 14,
      height      = 5,
      dpi         = 150,
    )
  } else if (plot_type == "circular") {
    CMplot(
      cm_df,
      plot.type   = "c",
      LOG10       = TRUE,
      threshold   = c(5e-8, 1e-5),
      file        = "png",
      file.name   = out_stem,
      file.output = TRUE,
      verbose     = FALSE,
      width       = 10,
      height      = 10,
      dpi         = 150,
    )
  } else if (plot_type == "qq") {
    CMplot(
      cm_df,
      plot.type   = "q",
      LOG10       = TRUE,
      file        = "png",
      file.name   = out_stem,
      file.output = TRUE,
      verbose     = FALSE,
      dpi         = 150,
    )
  }

  # Identify the file CMplot just wrote (new PNG not present before the call)
  after   <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
  created <- setdiff(after, before)

  if (length(created) == 0) {
    stop(sprintf("CMplot produced no PNG in %s (stem: %s)", out_dir, out_stem))
  }
  if (length(created) > 1) {
    warning(sprintf("CMplot wrote multiple PNGs; renaming first: %s",
                    paste(basename(created), collapse=", ")))
  }

  # Rename to the standardised path expected by the benchmark harness
  if (created[1] != out_path) {
    file.rename(created[1], out_path)
  }
}


run_qqman <- function(df, out_path, plot_type) {
  suppressPackageStartupMessages(library(qqman))

  png(out_path, width=1400, height=500, res=150)
  tryCatch({
    if (plot_type %in% c("manhattan", "circular")) {
      # qqman only does linear Manhattan; circular is not supported
      manhattan(
        df,
        chr  = "CHR",
        bp   = "BP",
        p    = "P",
        snp  = "SNP",
        main = ""
      )
    } else if (plot_type == "qq") {
      qq(df$P, main="")
    }
  }, finally = {
    dev.off()
  })
}



# ---------------------------------------------------------------------------
# run_CMplot_multitrack: multi-trait / multi-track CMplot benchmark
#
# CMplot requires traits to be merged column-wise into a single data frame
# before plotting — the merge is performed HERE as part of the timed section
# because it is required user work that has no equivalent in pycmplot.
# pycmplot simply accepts a list of file paths (see bench_python.py).
#
# Input format comparison (key finding for manuscript Table 1):
#   pycmplot : list of N separate TSV files — no preparation needed
#   CMplot   : one merged df with cols SNP, CHR, BP, Trait1, Trait2, ... TraitN
# ---------------------------------------------------------------------------
run_CMplot_multitrack <- function(manifest_path, out_path, plot_type) {
  suppressPackageStartupMessages(library(CMplot))

  # --- Read manifest and load trait files (timed) --------------------------
  traits_path <- readLines(manifest_path)
  trait_paths <- readLines(manifest_path)
  trait_paths <- trait_paths[nzchar(trimws(trait_paths))]
  n_traits    <- length(trait_paths)
  trait_labels <- paste0("Trait", seq_len(n_traits))

  # Load first file to get scaffold; load remaining p-values only
  base_df <- fread(trait_paths[1], sep="\t", select=c("SNP","CHR","BP","P"),
                   data.table=FALSE)
  colnames(base_df)[4] <- trait_labels[1]

  for (k in seq(2, n_traits)) {
    pvals <- fread(trait_paths[k], sep="\t", select="P", data.table=FALSE)[[1]]
    base_df[[trait_labels[k]]] <- pvals
  }
  # base_df now has: SNP, CHR, BP, Trait1, Trait2, ..., TraitN
  # This is the only input format CMplot accepts for multi-trait plots.
  
  #print(head(base_df))

  # --- CMplot call (timed — merge above is included) -----------------------
  out_dir  <- normalizePath(dirname(out_path), mustWork=TRUE)
  out_stem <- tools::file_path_sans_ext(basename(out_path))
  before   <- list.files(out_dir, pattern="\\.png$", full.names=TRUE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add=TRUE)
  setwd(out_dir)

  cm_type <- if (grepl("circular", plot_type)) "c" else "m"

  if (grepl("circular", plot_type)) {
    CMplot(
      base_df,
      plot.type   = 'c',
      LOG10       = TRUE,
      threshold   = c(5e-8, 1e-5),
      file        = "png",
      file.name   = out_stem,
      file.output = TRUE,
      verbose     = FALSE,
      width       = if (cm_type == "m") 14 else 10,
      height      = if (cm_type == "m") 5  else 10,
      dpi         = 150,
    )
  } else {
    CMplot(
      base_df,
      plot.type   = 'm',
      LOG10       = TRUE,
      threshold   = c(5e-8, 1e-5),
      file        = "png",
      file.name   = out_stem,
      file.output = TRUE,
      verbose     = FALSE,
      width       = if (cm_type == "m") 14 else 10,
      height      = if (cm_type == "m") 5  else 10,
      dpi         = 150,
      multracks   = TRUE,
    )
  }

  after   <- list.files(out_dir, pattern="\\.png$", full.names=TRUE)
  created <- setdiff(after, before)
  if (length(created) == 0) stop(sprintf("CMplot wrote no PNG (stem: %s)", out_stem))
  if (created[1] != out_path) file.rename(created[1], out_path)
}

RUNNERS <- list(
  CMplot              = run_CMplot,
  CMplot_multitrack   = run_CMplot_multitrack,
  qqman               = run_qqman
)


# ---------------------------------------------------------------------------
# Memory helper — peak RSS via gc() + object.size heuristics
# R has no clean tracemalloc equivalent; we use gc() before/after and
# report the difference in Vcells (heap cells, 8 bytes each on 64-bit).
# For a more accurate measure, the caller can wrap with /usr/bin/time -v.
# ---------------------------------------------------------------------------
peak_mem_mb <- function(expr) {
  gc(reset=TRUE, full=TRUE)
  mem_before <- sum(gc(reset=FALSE)[,2])   # Vcells in use
  force(expr)
  gc_after <- gc(reset=FALSE)
  mem_after <- sum(gc_after[,2])
  # Vcells are 8 bytes each
  delta_mb <- (mem_after - mem_before) * 8 / 1024 / 1024
  max(delta_mb, 0)
}


# ---------------------------------------------------------------------------
# Main benchmark loop
# ---------------------------------------------------------------------------

#bench_res <- paste0(tool, "_bench_r.csv")
bench_res <- paste0(tool, "_bench_r.csv")
csv_path <- file.path(outdir, bench_res)

cols <- c("tool","plot_type","size_label","n_variants",
          "replicate","wall_time_s","peak_mem_mb","out_file_kb")

# For multi-track tools, input_for_bench is the manifest path;
# for single-trait tools it is the TSV file path (--input).
is_multitrack <- !is.null(manifest)
input_for_bench <- if (is_multitrack) manifest else input

cat(sprintf("\n[bench] tool=%s  size=%s  plot=%s  reps=%d\n",
            tool, size, plot_type, reps))

# Load data for single-trait tools; for multitrack tools each runner handles
# its own I/O inside the timed section (merge is part of CMplot's cost).
if (!is_multitrack) {
  cat("  Loading data... ")
  df_full <- fread(input_for_bench, sep="\t", data.table=FALSE)
  n_variants <- nrow(df_full)
  cat(sprintf("%s variants\n", format(n_variants, big.mark=",")))
} else {
  # Report variant count from the first trait file in the manifest
  first_path <- readLines(manifest, n=1)
  n_variants <- nrow(fread(first_path, select="SNP", data.table=FALSE))
  n_traits   <- length(readLines(manifest)[nzchar(trimws(readLines(manifest)))])
  cat(sprintf("  Multi-track: %d traits x %s variants per trait\n",
              n_traits, format(n_variants, big.mark=",")))
  df_full <- NULL  # not used; runners load their own data
}

runner <- RUNNERS[[tool]]
if (is.null(runner)) stop(sprintf("Unknown tool: %s", tool))

# Open CSV (append mode)
write_header <- !file.exists(csv_path)
fh <- file(csv_path, open="a")
if (write_header) {
  writeLines(paste(cols, collapse=","), fh)
}

# When --replicate N is supplied (called from bench_r_timed.sh), run only
# that single replicate. Memory is written as EXTERNAL — the wrapper script
# captures /usr/bin/time -v MaxRSS and patches it into the CSV afterwards.
rep_seq <- if (!is.null(single_rep)) single_rep else seq_len(reps)
external_mem <- !is.null(single_rep)

for (rep in rep_seq) {
  out_fig <- file.path(figdir,
    sprintf("%s_%s_%s_rep%d.png", tool, plot_type, size, rep))

  # Fresh copy each replicate to avoid caching effects
  # For multitrack runs df_full is NULL; runners handle their own I/O.
  df <- if (!is_multitrack) df_full else NULL

  result <- tryCatch({
    gc(full=TRUE, reset=TRUE)
    mem_vcells_before <- sum(gc(reset=FALSE)[,2])

    t0 <- proc.time()[["elapsed"]]
    if (is_multitrack) {
      # Multi-track runners receive the manifest path; they handle all I/O
      # (including the column-wise merge for CMplot) inside the timed section.
      runner(input_for_bench, out_fig, plot_type)
    } else {
      runner(df, out_fig, plot_type)
    }
    t1 <- proc.time()[["elapsed"]]

    mem_vcells_after <- sum(gc(reset=FALSE)[,2])
    wall_s  <- round(t1 - t0, 3)

    # Use gc()-based estimate when running standalone; EXTERNAL when wrapped
    # by bench_r_timed.sh so /usr/bin/time -v MaxRSS can be patched in.
    mem_mb  <- if (external_mem) "EXTERNAL" else
               round(max((mem_vcells_after - mem_vcells_before) * 8 / 1024^2, 0), 2)
    file_kb <- if (file.exists(out_fig)) round(file.info(out_fig)$size / 1024, 1) else 0

    cat(sprintf("  rep %d  time=%.2fs  mem=%s  fig=%.0fKB\n",
                rep, wall_s,
                if (external_mem) "EXTERNAL" else sprintf("%.0fMB", as.numeric(mem_mb)),
                file_kb))

    c(tool, plot_type, size, n_variants, rep, wall_s, mem_mb, file_kb)
  }, error = function(e) {
    cat(sprintf("  rep %d FAILED: %s\n", rep, conditionMessage(e)), file=stderr())
    c(tool, plot_type, size, n_variants, rep, "ERROR", "ERROR", "ERROR")
  })

  writeLines(paste(result, collapse=","), fh)
  flush(fh)
}

close(fh)
cat(sprintf("\nResults appended to %s\n", csv_path))



# #!/usr/bin/env Rscript
# # bench_r.R
# # Benchmarks R GWAS visualization tools: CMplot, qqman.
# #
# # Usage:
# #   Rscript bench_r.R --tool CMplot --input data/sumstats_1M.tsv \
# #       --size 1M --replicates 5 --outdir results/ --figdir figures/
# #
# # Writes one CSV row per replicate to results/bench_r.csv.
# 
# suppressPackageStartupMessages({
#   library(optparse)
#   library(data.table)
# })
# 
# # ---------------------------------------------------------------------------
# # CLI arguments
# # ---------------------------------------------------------------------------
# option_list <- list(
#   make_option("--tool",       type="character", help="Tool: CMplot or qqman"),
#   make_option("--input",      type="character", help="Path to sumstats TSV"),
#   make_option("--size",       type="character", help="Dataset size label (e.g. 1M)"),
#   make_option("--plot-type",  type="character", default="manhattan",
#               help="Plot type: manhattan, circular, qq [default: manhattan]"),
#   make_option("--replicates", type="integer",   default=5L),
#   make_option("--replicate",  type="integer",   default=NULL,
#               help="Run only this single replicate (1-based). Used by bench_r_timed.sh\n\t\t    to isolate each rep in its own process for /usr/bin/time -v memory capture.\n\t\t    When set, peak_mem_mb is written as EXTERNAL (filled in by the wrapper)."),
#   make_option("--outdir",     type="character", default="results"),
#   make_option("--figdir",     type="character", default="figures")
# )
# 
# opt <- parse_args(OptionParser(option_list=option_list))
# tool       <- opt$tool
# input      <- opt$input
# size       <- opt$size
# plot_type  <- opt[["plot-type"]]
# reps       <- opt$replicates
# single_rep <- opt$replicate   # NULL unless --replicate N supplied
# outdir     <- opt$outdir
# figdir     <- opt$figdir
# 
# dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
# dir.create(figdir, showWarnings=FALSE, recursive=TRUE)
# 
# # ---------------------------------------------------------------------------
# # Tool wrappers
# # ---------------------------------------------------------------------------
# 
# run_CMplot <- function(df, out_path, plot_type) {
#   suppressPackageStartupMessages(library(CMplot))
# 
#   # CMplot expects: SNP, CHR, BP, then one or more p-value columns.
#   # This specific column order is a known usability friction point.
#   cm_df <- df[, c("SNP", "CHR", "BP", "P")]
# 
#   # --- CMplot output file naming fix ---------------------------------------
#   # CMplot constructs its own filename by prepending a plot-type prefix
#   # (e.g. "Rect_Manhtn.", "Circ_Manhtn.", "QQplot.") to whatever you pass
#   # as file.name, and writes the result relative to the *working directory*.
#   # Passing a full absolute path as file.name therefore produces an invalid
#   # concatenation like "Rect_Manhtn./full/path/stem" which cannot be opened.
#   #
#   # Fix: (1) cd into the figures directory, (2) pass only the bare stem,
#   # (3) use list.files() with a pattern to find whatever CMplot actually
#   #     wrote (prefix varies across CMplot versions), then rename to out_path.
#   # -------------------------------------------------------------------------
#   out_dir  <- normalizePath(dirname(out_path), mustWork = TRUE)
#   out_stem <- tools::file_path_sans_ext(basename(out_path))
# 
#   # Snapshot files present before calling CMplot so we can identify new ones
#   before <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
# 
#   # Temporarily change working directory; restore on function exit
#   old_wd <- getwd()
#   on.exit(setwd(old_wd), add = TRUE)
#   setwd(out_dir)
# 
#   if (plot_type == "manhattan") {
#     CMplot(
#       cm_df,
#       plot.type   = "m",
#       LOG10       = TRUE,
#       threshold   = c(5e-8, 1e-5),
#       file        = "png",
#       file.name   = out_stem,       # bare stem only — no path, no extension
#       file.output = TRUE,
#       verbose     = FALSE,
#       width       = 14,
#       height      = 5,
#       dpi         = 150,
#     )
#   } else if (plot_type == "circular") {
#     CMplot(
#       cm_df,
#       plot.type   = "c",
#       LOG10       = TRUE,
#       threshold   = c(5e-8, 1e-5),
#       file        = "png",
#       file.name   = out_stem,
#       file.output = TRUE,
#       verbose     = FALSE,
#       width       = 10,
#       height      = 10,
#       dpi         = 150,
#     )
#   } else if (plot_type == "qq") {
#     CMplot(
#       cm_df,
#       plot.type   = "q",
#       LOG10       = TRUE,
#       file        = "png",
#       file.name   = out_stem,
#       file.output = TRUE,
#       verbose     = FALSE,
#       dpi         = 150,
#     )
#   }
# 
#   # Identify the file CMplot just wrote (new PNG not present before the call)
#   after   <- list.files(out_dir, pattern = "\\.png$", full.names = TRUE)
#   created <- setdiff(after, before)
# 
#   if (length(created) == 0) {
#     stop(sprintf("CMplot produced no PNG in %s (stem: %s)", out_dir, out_stem))
#   }
#   if (length(created) > 1) {
#     warning(sprintf("CMplot wrote multiple PNGs; renaming first: %s",
#                     paste(basename(created), collapse=", ")))
#   }
# 
#   # Rename to the standardised path expected by the benchmark harness
#   if (created[1] != out_path) {
#     file.rename(created[1], out_path)
#   }
# }
# 
# 
# run_qqman <- function(df, out_path, plot_type) {
#   suppressPackageStartupMessages(library(qqman))
# 
#   png(out_path, width=1400, height=500, res=150)
#   tryCatch({
#     if (plot_type %in% c("manhattan", "circular")) {
#       # qqman only does linear Manhattan; circular is not supported
#       manhattan(
#         df,
#         chr  = "CHR",
#         bp   = "BP",
#         p    = "P",
#         snp  = "SNP",
#         main = ""
#       )
#     } else if (plot_type == "qq") {
#       qq(df$P, main="")
#     }
#   }, finally = {
#     dev.off()
#   })
# }
# 
# 
# RUNNERS <- list(
#   CMplot = run_CMplot,
#   qqman  = run_qqman
# )
# 
# 
# # ---------------------------------------------------------------------------
# # Memory helper — peak RSS via gc() + object.size heuristics
# # R has no clean tracemalloc equivalent; we use gc() before/after and
# # report the difference in Vcells (heap cells, 8 bytes each on 64-bit).
# # For a more accurate measure, the caller can wrap with /usr/bin/time -v.
# # ---------------------------------------------------------------------------
# peak_mem_mb <- function(expr) {
#   gc(reset=TRUE, full=TRUE)
#   mem_before <- sum(gc(reset=FALSE)[,2])   # Vcells in use
#   force(expr)
#   gc_after <- gc(reset=FALSE)
#   mem_after <- sum(gc_after[,2])
#   # Vcells are 8 bytes each
#   delta_mb <- (mem_after - mem_before) * 8 / 1024 / 1024
#   max(delta_mb, 0)
# }
# 
# 
# # ---------------------------------------------------------------------------
# # Main benchmark loop
# # ---------------------------------------------------------------------------
# 
# bench_res <- paste0(tool, "_bench_r.csv")
# 
# csv_path <- file.path(outdir, bench_res)
# cols <- c("tool","plot_type","size_label","n_variants",
#           "replicate","wall_time_s","peak_mem_mb","out_file_kb")
# 
# cat(sprintf("\n[bench] tool=%s  size=%s  plot=%s  reps=%d\n",
#             tool, size, plot_type, reps))
# 
# # Load data once outside the replicate loop — we time each plot call
# # including a fresh subset copy so memory pressure is realistic.
# cat("  Loading data... ")
# df_full <- fread(input, sep="\t", data.table=FALSE)
# n_variants <- nrow(df_full)
# cat(sprintf("%s variants\n", format(n_variants, big.mark=",")))
# 
# runner <- RUNNERS[[tool]]
# if (is.null(runner)) stop(sprintf("Unknown tool: %s", tool))
# 
# # Open CSV (append mode)
# write_header <- !file.exists(csv_path)
# fh <- file(csv_path, open="a")
# if (write_header) {
#   writeLines(paste(cols, collapse=","), fh)
# }
# 
# # When --replicate N is supplied (called from bench_r_timed.sh), run only
# # that single replicate. Memory is written as EXTERNAL — the wrapper script
# # captures /usr/bin/time -v MaxRSS and patches it into the CSV afterwards.
# rep_seq <- if (!is.null(single_rep)) single_rep else seq_len(reps)
# external_mem <- !is.null(single_rep)
# 
# for (rep in rep_seq) {
#   out_fig <- file.path(figdir,
#     sprintf("%s_%s_%s_rep%d.png", tool, plot_type, size, rep))
# 
#   # Fresh copy each replicate to avoid caching effects
#   df <- df_full
# 
#   result <- tryCatch({
#     gc(full=TRUE, reset=TRUE)
#     mem_vcells_before <- sum(gc(reset=FALSE)[,2])
# 
#     t0 <- proc.time()[["elapsed"]]
#     runner(df, out_fig, plot_type)
#     t1 <- proc.time()[["elapsed"]]
# 
#     mem_vcells_after <- sum(gc(reset=FALSE)[,2])
#     wall_s  <- round(t1 - t0, 3)
# 
#     # Use gc()-based estimate when running standalone; EXTERNAL when wrapped
#     # by bench_r_timed.sh so /usr/bin/time -v MaxRSS can be patched in.
#     mem_mb  <- if (external_mem) "EXTERNAL" else
#                round(max((mem_vcells_after - mem_vcells_before) * 8 / 1024^2, 0), 2)
#     file_kb <- if (file.exists(out_fig)) round(file.info(out_fig)$size / 1024, 1) else 0
# 
#     cat(sprintf("  rep %d  time=%.2fs  mem=%s  fig=%.0fKB\n",
#                 rep, wall_s,
#                 if (external_mem) "EXTERNAL" else sprintf("%.0fMB", as.numeric(mem_mb)),
#                 file_kb))
# 
#     c(tool, plot_type, size, n_variants, rep, wall_s, mem_mb, file_kb)
#   }, error = function(e) {
#     cat(sprintf("  rep %d FAILED: %s\n", rep, conditionMessage(e)), file=stderr())
#     c(tool, plot_type, size, n_variants, rep, "ERROR", "ERROR", "ERROR")
#   })
# 
#   writeLines(paste(result, collapse=","), fh)
#   flush(fh)
# }
# 
# close(fh)
# cat(sprintf("\nResults appended to %s\n", csv_path))

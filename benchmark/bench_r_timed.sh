#!/usr/bin/env bash
# bench_r_timed.sh
#
# Runs each R benchmark replicate in its own Rscript process wrapped with
# /usr/bin/time -v, capturing OS-level peak RSS (Maximum resident set size)
# per replicate. Replaces the gc()-based memory estimate in bench_r.R with
# the more accurate /usr/bin/time -v value.
#
# Usage (mirrors bench_r.R arguments):
#   bash bench_r_timed.sh \
#       --tool CMplot --input data/sumstats_1M.tsv \
#       --size 1M --plot-type manhattan \
#       --replicates 5 --outdir results/ --figdir figures/
#
# Requires:
#   /usr/bin/time (GNU time, not bash built-in — ships with util-linux on Linux)
#   Rscript with CMplot / qqman installed
#
# Output:
#   Appends to results/bench_r.csv (same file as bench_r.R standalone mode)
#   with peak_mem_mb filled from /usr/bin/time -v MaxRSS rather than gc().
# --------------------------------------------------------------------------

set -euo pipefail

# --------------------------------------------------------------------------
# Defaults (overridden by flags below)
# --------------------------------------------------------------------------
TOOL=""
INPUT=""
SIZE=""
MANIFEST=""
PLOT_TYPE="manhattan"
REPLICATES=5
OUTDIR="results"
FIGDIR="figures"

# --------------------------------------------------------------------------
# Parse arguments (same interface as bench_r.R)
# --------------------------------------------------------------------------
while [[ $# -gt 0 ]]; do
  case "$1" in
    --tool)        TOOL="$2";       shift 2 ;;
    --input)       INPUT="$2";      shift 2 ;;
    --size)        SIZE="$2";       shift 2 ;;
    --plot-type)   PLOT_TYPE="$2";  shift 2 ;;
    --replicates)  REPLICATES="$2"; shift 2 ;;
    --manifest)    MANIFEST="$2";   shift 2 ;;
    --outdir)      OUTDIR="$2";     shift 2 ;;
    --figdir)      FIGDIR="$2";     shift 2 ;;
    *) echo "Unknown argument: $1" >&2; exit 1 ;;
  esac
done

#if [[ -z "$TOOL" || -z "$INPUT" || -z "$SIZE" ]]; then
#  echo "Usage: bash bench_r_timed.sh --tool <tool> --input <file> --size <label> [options]" >&2
#  exit 1
#fi

if [[ -z "$TOOL" || -z "$SIZE" ]]; then
  echo "Usage: bash bench_r_timed.sh --tool <tool> --size <label> [options]" >&2
  exit 1
fi

# --------------------------------------------------------------------------
# Verify /usr/bin/time is GNU time (not the bash built-in)
# --------------------------------------------------------------------------
if ! /usr/bin/time --version 2>&1 | grep -q "GNU"; then
  echo "ERROR: /usr/bin/time is not GNU time." >&2
  echo "  On macOS install via: brew install gnu-time" >&2
  echo "  then use gtime instead of /usr/bin/time." >&2
  exit 1
fi

# --------------------------------------------------------------------------
# Directory setup
# --------------------------------------------------------------------------
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
LOGDIR="${WORKDIR}/logs"
mkdir -p "$OUTDIR" "$FIGDIR" "$LOGDIR"

# --------------------------------------------------------------------------
# CSV setup — same columns as bench_r.R
# --------------------------------------------------------------------------
CSV="${OUTDIR}/${TOOL}_bench_r.csv"
COLS="tool,plot_type,size_label,n_variants,replicate,wall_time_s,peak_mem_mb,out_file_kb"
[[ ! -f "$CSV" ]] && echo "$COLS" > "$CSV"

# --------------------------------------------------------------------------
# parse_maxrss: extract Maximum resident set size (kbytes) from
# /usr/bin/time -v stderr, convert to MB
# --------------------------------------------------------------------------
parse_maxrss() {
  local time_log="$1"
  local kb
  kb=$(grep "Maximum resident set size" "$time_log" | awk '{print $NF}')
  if [[ -z "$kb" || "$kb" == "0" ]]; then
    echo "NA"
  else
    # Convert kbytes -> MB, two decimal places
    awk "BEGIN { printf \"%.2f\", ${kb} / 1024 }"
  fi
}

# --------------------------------------------------------------------------
# patch_external: replace the EXTERNAL placeholder written by bench_r.R
# with the actual MaxRSS value. Targets the last occurrence of EXTERNAL
# in the CSV (the row just written by this replicate's Rscript call).
# Uses a temp file to avoid in-place sed portability issues.
# --------------------------------------------------------------------------
patch_external() {
  local csv="$1"
  local maxrss_mb="$2"
  local tmpfile
  tmpfile="$(mktemp)"
  # Reverse the file, replace first EXTERNAL (= last in file), reverse back
  awk -v val="$maxrss_mb" '
    !done && /EXTERNAL/ { sub(/EXTERNAL/, val); done=1 }
    { print }
  ' <(tac "$csv") | tac > "$tmpfile"
  mv "$tmpfile" "$csv"
}

# --------------------------------------------------------------------------
# Main loop — one Rscript process per replicate
# --------------------------------------------------------------------------
echo ""
echo "[timed] tool=${TOOL}  size=${SIZE}  plot=${PLOT_TYPE}  reps=${REPLICATES}"
echo "[timed] Memory source: /usr/bin/time -v MaxRSS"
echo ""

for REP in $(seq 1 "$REPLICATES"); do
  TIME_LOG="${LOGDIR}/${TOOL}_${PLOT_TYPE}_${SIZE}_rep${REP}.timev"

  echo -n "  rep ${REP}/${REPLICATES} ... "

  # Run bench_r.r for this single replicate, capturing /usr/bin/time -v
  # stderr to TIME_LOG. bench_r.R stdout goes to the terminal as normal.
  # Note: /usr/bin/time -v writes its report to stderr; we redirect only
  # that to TIME_LOG while letting Rscript's own stderr (warnings etc.)
  # pass through to the terminal.
  #set +e
  #/usr/bin/time -v \
  #  Rscript "${WORKDIR}/bench_r.r" \
  #    --tool       "$TOOL" \
  #    --input      "$INPUT" \
  #    --size       "$SIZE" \
  #    --manifest   "$MANIFEST" \
  #    --plot-type  "$PLOT_TYPE" \
  #    --replicate  "$REP" \
  #    --outdir     "$OUTDIR" \
  #    --figdir     "$FIGDIR" \
  #  2> "$TIME_LOG"
  #EXIT_CODE=$?
  #set -e

  set +e
  /usr/bin/time -v \
    Rscript "${WORKDIR}/bench_r.r" \
      --tool       "$TOOL" \
      --size       "$SIZE" \
      --manifest   "$MANIFEST" \
      --plot-type  "$PLOT_TYPE" \
      --replicate  "$REP" \
      --outdir     "$OUTDIR" \
      --figdir     "$FIGDIR" \
    2> "$TIME_LOG"
  EXIT_CODE=$?
  set -e

  if [[ $EXIT_CODE -ne 0 ]]; then
    echo "FAILED (exit ${EXIT_CODE})"
    # bench_r.r writes its own ERROR row; nothing to patch
    continue
  fi

  # Parse MaxRSS and patch the EXTERNAL placeholder in the CSV
  MAXRSS_MB=$(parse_maxrss "$TIME_LOG")
  patch_external "$CSV" "$MAXRSS_MB"

  echo "MaxRSS=${MAXRSS_MB}MB  (log: $(basename "$TIME_LOG"))"
done

echo ""
echo "[timed] Done. Results in: ${CSV}"

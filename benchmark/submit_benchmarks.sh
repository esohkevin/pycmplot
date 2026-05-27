#!/usr/bin/env bash
# submit_benchmarks.sh
#
# Submits one SLURM job per tool, each running all dataset sizes and
# replicates sequentially on an identical, isolated compute node.
#
# Usage:
#   bash submit_benchmarks.sh
#   bash submit_benchmarks.sh --dry-run      # print commands without submitting
#   bash submit_benchmarks.sh --sizes "500K 1M"  # override dataset sizes
#
# Requirements:
#   - SLURM environment
#   - Conda environments: gwas_bench_py (Python tools), gwas_bench_r (R tools)
#   - Data directory populated by: python generate_sumstats.py --all
# -----------------------------------------------------------------------

set -euo pipefail

# -----------------------------------------------------------------------
# Configuration — edit these for your cluster
# -----------------------------------------------------------------------
WORKDIR="$(cd "$(dirname "$0")" && pwd)"
DATADIR="${WORKDIR}/data"
DATADIR="/scratch4/awonkam1/kesoh/gwas/pycmplotbenchmark/"
RESULTSDIR="${WORKDIR}/results"
FIGDIR="${WORKDIR}/figures"
LOGDIR="${WORKDIR}/logs"

PARTITION="shared"           # your cluster's partition name
TIME_LIMIT="12:00:00"      # generous limit; large R runs can be slow
MEM_PER_NODE="64G"         # enough for 10M-variant datasets
CPUS=1                     # SINGLE CORE — critical for fair comparison
CONDA_PY="gwas_bench_py"   # conda env with pycmplot, gwaslab, qmplot
CONDA_R="gwas_bench_r"     # conda env with R + CMplot, qqman
REPLICATES=5

# Dataset sizes to benchmark (must match files in DATADIR)
SIZES=("500K" "1M" "2M" "5M" "10M")

# Plot types per tool (space-separated)
declare -A PLOT_TYPES
PLOT_TYPES[pycmplot]="manhattan circular qq"
PLOT_TYPES[gwaslab]="manhattan qq"
PLOT_TYPES[qmplot]="manhattan qq"
PLOT_TYPES[CMplot]="manhattan circular qq"
PLOT_TYPES[qqman]="manhattan qq"

# -----------------------------------------------------------------------
# Parse flags
# -----------------------------------------------------------------------
DRY_RUN=0
for arg in "$@"; do
  case "$arg" in
    --dry-run) DRY_RUN=1 ;;
    --sizes=*) IFS=' ' read -ra SIZES <<< "${arg#--sizes=}" ;;
  esac
done

mkdir -p "$RESULTSDIR" "$FIGDIR" "$LOGDIR"

# -----------------------------------------------------------------------
# Helper: submit one SLURM job
# -----------------------------------------------------------------------
submit_job() {
  local job_name="$1"
  local script_body="$2"

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "--- [DRY RUN] $job_name ---"
    echo "$script_body"
    echo ""
    return
  fi

  local job_id
  job_id=$(echo "$script_body" | sbatch --parsable)
  echo "Submitted $job_name  -> job $job_id"
}

# -----------------------------------------------------------------------
# Python tools
# -----------------------------------------------------------------------
for TOOL in pycmplot gwaslab qmplot; do

  SCRIPT=$(cat <<SLURM
#!/usr/bin/env bash
#SBATCH --job-name=bench_${TOOL}
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH -w sr47
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM_PER_NODE}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGDIR}/bench_${TOOL}_%j.out
#SBATCH --error=${LOGDIR}/bench_${TOOL}_%j.err

# -- Environment isolation: single core, no implicit threading
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

# -- Activate Python environment
#source "\$(conda info --base)/etc/profile.d/conda.sh"
module load anaconda
conda activate ${CONDA_PY}

echo "=============================="
echo "Job: bench_${TOOL}"
echo "Node: \$(hostname)"
echo "Date: \$(date)"
echo "Python: \$(python --version)"
echo "=============================="

cd ${WORKDIR}

for SIZE in ${SIZES[@]}; do
  INPUT="${DATADIR}/sumstats_\${SIZE}.tsv"
  if [[ ! -f "\$INPUT" ]]; then
    echo "WARNING: \$INPUT not found — skipping"
    continue
  fi
  for PLOT_TYPE in ${PLOT_TYPES[$TOOL]}; do
    echo ""
    echo ">>> ${TOOL}  size=\${SIZE}  plot=\${PLOT_TYPE}"
    python bench_python.py \
      --tool ${TOOL} \
      --input "\$INPUT" \
      --size "\${SIZE}" \
      --plot-type "\${PLOT_TYPE}" \
      --replicates ${REPLICATES} \
      --outdir "${RESULTSDIR}" \
      --figdir "${FIGDIR}"
  done
done

echo ""
echo "Done: \$(date)"
SLURM
)

  submit_job "bench_${TOOL}" "$SCRIPT"
done

# -----------------------------------------------------------------------
# R tools
# -----------------------------------------------------------------------
for TOOL in CMplot qqman; do

  SCRIPT=$(cat <<SLURM
#!/usr/bin/env bash
#SBATCH --job-name=bench_${TOOL}
#SBATCH --partition=${PARTITION}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=${CPUS}
#SBATCH --mem=${MEM_PER_NODE}
#SBATCH --time=${TIME_LIMIT}
#SBATCH --output=${LOGDIR}/bench_${TOOL}_%j.out
#SBATCH --error=${LOGDIR}/bench_${TOOL}_%j.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

source "\$(conda info --base)/etc/profile.d/conda.sh"
conda activate ${CONDA_R}

echo "=============================="
echo "Job: bench_${TOOL}"
echo "Node: \$(hostname)"
echo "Date: \$(date)"
echo "R: \$(R --version | head -1)"
echo "=============================="

cd ${WORKDIR}

for SIZE in ${SIZES[@]}; do
  INPUT="${DATADIR}/sumstats_\${SIZE}.tsv"
  if [[ ! -f "\$INPUT" ]]; then
    echo "WARNING: \$INPUT not found — skipping"
    continue
  fi
  for PLOT_TYPE in ${PLOT_TYPES[$TOOL]}; do
    echo ""
    echo ">>> ${TOOL}  size=\${SIZE}  plot=\${PLOT_TYPE}"
    Rscript bench_r.R \
      --tool ${TOOL} \
      --input "\$INPUT" \
      --size "\${SIZE}" \
      --plot-type "\${PLOT_TYPE}" \
      --replicates ${REPLICATES} \
      --outdir "${RESULTSDIR}" \
      --figdir "${FIGDIR}"
  done
done

echo ""
echo "Done: \$(date)"
SLURM
)

  submit_job "bench_${TOOL}" "$SCRIPT"
done

echo ""
echo "All jobs submitted. Monitor with: squeue -u \$USER"
echo "Results will accumulate in: ${RESULTSDIR}/"

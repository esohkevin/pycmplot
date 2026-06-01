#!/usr/bin/env bash
#SBATCH --job-name=bench_qqman
#SBATCH --partition=shared
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH -w sr47
#SBATCH --time=12:00:00
#SBATCH --output=logs/bench_qqman_%j.out
#SBATCH --error=logs/bench_qqman_%j.err

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export R_DATATABLE_NUM_THREADS=1

#source "$(conda info --base)/etc/profile.d/conda.sh"
module load anaconda
conda activate gwas_bench_r

mkdir -p logs

echo "=============================="
echo "Job: bench_qqman"
echo "Node: $(hostname)"
echo "Date: $(date)"
echo "R: $(R --version | head -1)"
echo "=============================="


OUTDIR="results"
FIGDIR="figures"
REPLICATES=5

for SIZE in 500K 1M 2M 5M 10M; do
  INPUT="data/sumstats_${SIZE}.tsv"
  if [[ ! -f "$INPUT" ]]; then
    echo "WARNING: $INPUT not found <E2><80><94> skipping"
    continue
  fi
  for PLOT_TYPE in manhattan circular qq; do
    echo ""
    echo ">>> qqman  size=${SIZE}  plot=${PLOT_TYPE}"

    bash bench_r_timed.sh \
      --tool qqman --input "$INPUT" --size "$SIZE" \
      --plot-type "$PLOT_TYPE" --replicates "$REPLICATES" \
      --outdir "$OUTDIR" --figdir "$FIGDIR"
  done
done

echo ""
echo "Done: $(date)"
